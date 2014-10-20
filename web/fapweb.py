"""
Faps web interface.

Show all available faps options as a webpage to help make fap files and submit
jobs. Built on flask, so it will run a webserver and tell the user how
to access the site.
"""

import argparse
import re
import socket
import webbrowser
from collections import defaultdict
from os import path
from subprocess import Popen, PIPE

from flask import Flask, request, render_template, Response, url_for, redirect

FAPS_ROOT = path.dirname(path.dirname(path.realpath(__file__)))


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
                        help="Do not try to launch a browser on start")
    args = parser.parse_args()
    return args


def parse_defaults():
    """Read all the options from the defaults.ini. Return as a dict."""
    defaults = open(path.join(FAPS_ROOT, 'defaults.ini'))

    for _header in xrange(12):
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


# create the application
# We run under /faps so that the URL is already prefixed for proxy'd
# connections
app = Flask(__name__, static_url_path='/faps/static')
app.config.from_object(__name__)


# If not running under a proxy, can just divert straight to the UI
@app.route('/')
def redirect_to_faps():
    """Root URL redirects to the application sitting at '/faps"""
    return redirect(url_for('faps'), code=302)


@app.route('/faps')
def faps():
    """Return the main UI."""
    return render_template('options.html', options=parse_defaults())


@app.route('/faps/submit', methods=['POST', 'GET', 'PUT'])
def submit_job():
    """
    Take whatever comes in through the form and run the faps job.

    We use the faps.py in the parent directory and use `python` to run it.
    """

    # Put the cif on disk
    # TODO: check for overwriting
    cif_file = request.files['cif-file']
    if not cif_file.filename.endswith('.cif'):
        return 1
    cif_file.save(cif_file.filename)

    basename = cif_file.filename[:-4]
    fap_file = request.form['fap-file']
    with open("{}.fap".format(basename), 'w') as fap_out:
        fap_out.write(fap_file)

    faps_script = path.join(FAPS_ROOT, 'faps.py')
    faps_command = ['python', faps_script, basename]
    faps_run = Popen(faps_command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = faps_run.communicate()
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

    print(response.__dict__)
    return response


def main():
    """
    Launch the flask application and any user defined tasks. Put stuff here
    That should not happen when running behind a web server.

    """

    args = commandline()
    hostname = socket.getfqdn()
    # Pick a random portm and hope that is doesn't get taken
    # between closing it and starting the app
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(('localhost', 0))
    desired_port = sock.getsockname()[1]
    sock.close()

    url = "http://{}:{}/faps".format(hostname, desired_port)
    if not args.no_browser:
        print("Attempting to start a web browser")
        # Launch after 1 second
        webbrowser.open(url)
    # App runs here
    app.run(debug=True, host=hostname, port=desired_port)


if __name__ == '__main__':
    main()
