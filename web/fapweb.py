"""
Faps web interface.

Show all available options as a webpage to help make fap files and submit
jobs. Built on flask, so it will run a webserver and tell the user how
to access the site.
"""

import re
from collections import defaultdict
from os import path
from subprocess import Popen, PIPE

from flask import Flask, request, render_template, Response, url_for, redirect


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


def parse_defaults():
    """Read all the options from the defaults.ini. Return as a dict."""
    defaults = open('../defaults.ini')

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

    faps_script = path.join(path.dirname(path.dirname(path.realpath(__file__))), 'faps.py')
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


if __name__ == '__main__':
    app.run(debug=True, host="0.0.0.0")
