#!/usr/bin/env python

from __future__ import unicode_literals

"""
Faps web interface.

The faps web interface is a tool to show all available faps options as an
interactive webpage to help make fap files and submit jobs.

To start the application just run:

    python fapweb.py

You can also add the --help option to get more information on the script.

"""

import argparse
import re
import socket
import webbrowser
from collections import defaultdict
from io import StringIO
from itertools import count
from os import path
from subprocess import Popen, PIPE
# Python 2/3 fix
try:
    from configparser import SafeConfigParser
except ImportError:
    from ConfigParser import SafeConfigParser

from flask import Flask, request, render_template, Response, url_for, redirect

FAPS_ROOT = path.dirname(path.dirname(path.realpath(__file__)))
DOT_FAPS_DIR = path.join(path.expanduser('~'), '.faps')


class Option(object):
    """Structure to hold information on an individual option."""
    def __init__(self):
        """Put defaults for things like advanced."""
        self.option_name = ""
        self.default = ""
        self.help_text = ""
        self.option_type = ""
        self.advanced = False
        self.section = ""


def commandline():
    """Parse commandline arguments and return the arguments object."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-browser', '-n', action='store_true',
                        help="Do not try to launch a browser on start.")
    parser.add_argument('--port', '-p', type=int,
                        help="Run on a specific port.")
    parser.add_argument('--debug', '-d', action='store_true',
                        help="Run the app in debugging mode.")
    args = parser.parse_args()
    return args


def parse_site_defaults():
    """Find where the script is and load defaults"""
    site_ini_path = path.join(FAPS_ROOT, 'site.ini')
    try:
        # Add file contents to empty string to upcast to
        # unicode in Python 2.7 since io needs unicode input.
        filetemp = open(site_ini_path, 'r')
        site_ini = '' + filetemp.read()
        filetemp.close()
        if not '[site_config]' in site_ini.lower():
            site_ini = '[site_config]\n' + site_ini
        site_ini = StringIO(site_ini)
    except IOError:
        # file does not exist so we just use a blank string
        site_ini = StringIO('[site_config]\n')

    site_ini_config = SafeConfigParser()
    site_ini_config.readfp(site_ini)
    return site_ini_config


def parse_defaults():
    """Read all the options from the defaults.ini. Return as a dict."""
    defaults = open(path.join(FAPS_ROOT, 'defaults.ini'))
    site_ini = parse_site_defaults()

    for _header in range(12):
        defaults.readline()

    options = defaultdict(list)
    current_option = Option()

    for line in defaults:
        line = line.strip()
        if not line:
            # blank line after, ignore
            continue
        elif line[0] in ['#', ';']:
            line = line.lstrip('#; \t')
            if line.startswith('@type'):
                current_option.option_type = line.split(None, 1)[1]
            elif line.startswith('@section'):
                current_option.section = line.split(None, 1)[1]
            elif line.startswith('@advanced'):
                current_option.advanced = True
            else:
                current_option.help_text += " %s\n" % line.strip()
        else:
            # Names are simple
            option_name, default = line.split('=')
            current_option.option_name = option_name.strip()
            if site_ini.has_option('site_config', current_option.option_name):
                default = site_ini.get('site_config',
                                       current_option.option_name)

            # Deal with specific types, bools can be real bools
            if current_option.option_type == 'bool':
                if default.strip().lower() in ["1", "yes", "true", "on"]:
                    current_option.default = True
                elif default.strip().lower() in ["0", "no", "false", "off"]:
                    current_option.default = False
            elif 'enum' in current_option.option_type:
                enum_values = re.split(r'[\s,{}]*', current_option.option_type)
                enum_values = [x for x in enum_values if x]
                current_option.option_type = tuple(enum_values)
                current_option.default = default.strip()
            else:
                current_option.default = default.strip()

            # Separate sections for advanced options
            if current_option.advanced:
                advanced_name = "{}_advanced".format(current_option.section)
                options[advanced_name].append(current_option)
            else:
                options[current_option.section].append(current_option)

            # start fresh for next option
            current_option = Option()

    return options


def parse_guests():
    """Get the names and descriptions of all available guests."""
    guest_list = []

    guests_libs = [path.join(FAPS_ROOT, 'guests.lib'),
                   path.join(DOT_FAPS_DIR, 'guests.lib'),
                   path.join('.', 'guests.lib')]
    guests = SafeConfigParser()
    guests.read(guests_libs)

    for guest_id in guests.sections():
        guest_info = (guest_id,
                      guests.get(guest_id, 'name'),
                      guests.get(guest_id, 'source'))
        if guest_id != 'template':
            guest_list.append(guest_info)

    return guest_list


# create the application
# We run under /faps so that the URL is already prefixed for proxy'd
# connections
app = Flask(__name__, static_url_path='/faps/static',
            template_folder=path.join(FAPS_ROOT, 'web', 'templates'),
            static_folder=path.join(FAPS_ROOT, 'web', 'static'))
app.config.from_object(__name__)


# If not running under a proxy, can just divert straight to the UI
@app.route('/')
def redirect_to_faps():
    """Root URL redirects to the application sitting at '/faps"""
    return redirect(url_for('faps'), code=302)


@app.route('/faps')
def faps():
    """Return the main UI."""
    return render_template('options.html', options=parse_defaults(),
                           guests=parse_guests())


@app.route('/faps/submit', methods=['POST', 'GET', 'PUT'])
def submit_job():
    """
    Take whatever comes in through the form and run the faps job.

    We use the faps.py in the parent directory and use `python` to run it.
    """

    # Put the cif on disk
    cif_file = request.files['cif-file']
    cif_filename = cif_file.filename
    if not cif_filename.endswith('.cif'):
        return 1

    cif_basename = cif_filename[:-4]

    # Do not overwrite an existing file
    if path.exists(cif_filename):
        for idx in count():
            if not path.exists("{}_{}.cif".format(cif_basename, idx)):
                cif_basename = "{}_{}".format(cif_basename, idx)
                break

    cif_file.save("{}.cif".format(cif_basename))

    # Extract the fap file text and put it on disk
    fap_file = request.form['fap-file']
    with open("{}.fap".format(cif_basename), 'w') as fap_out:
        fap_out.write(fap_file)

    faps_script = path.join(FAPS_ROOT, 'faps.py')
    faps_command = ['python', faps_script, cif_basename]
    faps_run = Popen(faps_command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = faps_run.communicate()
    # in python3 these are byte arrays, need them as strings
    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')
    print('---err---')
    print(stderr)
    print('---out---')
    print(stdout)

    if "Faps terminated normally" in stdout:
        response = Response("Submitted successfully",
                            content_type='text/xml; charset=utf-8')
    else:
        response = Response("Something failed :(",
                            content_type='text/xml; charset=utf-8')

    print("Sending reply: {}".format(response.response))
    return response


def main():
    """
    Launch the flask application and any user defined tasks. Put stuff here
    That should not happen when running behind a web server.

    """

    args = commandline()

    # Information about how and where to run
    hostname = socket.getfqdn()
    if args.port:
        # If this is in use, just let everything fail.
        # If you've specified a port, then you should know what you are doing.
        desired_port = args.port
    else:
        # Pick a random port and hope that is doesn't get taken
        # between closing it and starting the app
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.bind(('localhost', 0))
        desired_port = sock.getsockname()[1]
        sock.close()

    url = "http://{}:{}/faps".format(hostname, desired_port)

    print(">> Starting faps web interface.")
    print(">> The UI is available at the URL: {}".format(url))
    print(">> The URL is accessible from any computer on the same network.")
    print(">> Use the browser on the machine with your cif files.")

    if args.no_browser:
        print(">> Web browser will not be launched. Go to the URL manually.")
    else:
        print(">> Attempting to start a web browser...")
        # Launch first and hope that the server is up by the time the
        # request is finished
        webbrowser.open(url)

    print(">> To exit press Ctrl-C")

    # App runs here...
    # Debugging mode uses the internal server, otherwise make use of
    # gevent to run things
    if args.debug:
        app.run(debug=True, host=hostname, port=desired_port)
    else:
        try:
            from gevent.wsgi import WSGIServer
            http_server = WSGIServer(listener=('', desired_port),
                                     application=app, log=None)
            http_server.serve_forever()
        except ImportError:
            from tornado.wsgi import WSGIContainer
            from tornado.httpserver import HTTPServer
            from tornado.ioloop import IOLoop
            from tornado.options import options, parse_command_line

            parse_command_line([None, '--logging=error'])

            http_server = HTTPServer(WSGIContainer(app))
            http_server.listen(desired_port)
            IOLoop.instance().start()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print(">> Exiting...")
        raise SystemExit
