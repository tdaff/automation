#! /usr/bin/python2.7
import glob
import os
from os import path
import json
import pickle
import random

import tornado.web
import tornado.ioloop
import tornado.template

from function_switch import ModifiableStructure, FunctionalGroupLibrary
from function_switch import random_combination_replace
from backend.web_store import WebStoreBackend


tornadoPort = 8888
cwd = os.getcwd() # used by static file server

DOT_FAPSWITCH_VERSION = (6, 0, 'w')


random_page = """
<!DOCTYPE html>
<html style="height:100%">
<title>HTML5-only and WebGL Jmol Demo</title>
<head>
<script type="text/javascript" src="js/JSmoljQuery.js"></script>
<script type="text/javascript" src="js/JSmolCore.js"></script>
<script type="text/javascript" src="js/JSmolApplet.js"></script>
<script type="text/javascript" src="js/JSmolApi.js"></script>
<script type="text/javascript" src="js/j2sjmol.js"></script>
<script type="text/javascript" src="js/JSmol.js"></script>
<!-- // following two only necessary for WebGL version: -->
<script type="text/javascript" src="js/JSmolThree.js"></script>
<script type="text/javascript" src="js/JSmolGLmol.js"></script>

  <script type="text/javascript">



// note that the variable name MUST match the first parameter in quotes


;(function() {{

// logic is set by indicating order of USE -- default is HTML5 for this test page, though
var use = "HTML5"//"JAVA HTML5 IMAGE" // WEBGL only by request by link below.
var s = document.location.search;


// Developers: The debugCode flag is checked in j2s/java/core.z.js,
// and, if TRUE, skips loading the core methods, forcing those
// to be read from their individual directories. Set this
// true if you want to do some code debugging by inserting
// System.out.println, document.title, or alert commands
// anywhere in the Java or Jmol code.

Jmol.debugCode = (s.indexOf("debugcode") >= 0);

//if (s.indexOf("?") < 0) s = "USE=HTML5"

if (s.indexOf("USE=") >= 0)
  use = s.split("USE=")[1].split("&")[0]
else if (s.indexOf("JAVA") >= 0)
  use = "JAVA"
else if (s.indexOf("IMAGE") >= 0)
  use = "IMAGE"
else if (s.indexOf("NOWEBGL") >= 0)
  use = "JAVA IMAGE"
else if (s.indexOf("WEBGL") >= 0)
  use = "WEBGL HTML5"
if (s.indexOf("NOWEBGL") >= 0)
  use = use.replace(/WEBGL/,"")
var useSignedApplet = (s.indexOf("SIGNED") >= 0);
if(useSignedApplet && use == "HTML5") use = "JAVA";

jmol_isReady = function(applet) {{
	document.title = (applet._id + " is ready")
	Jmol._getElement(applet, "appletdiv").style.border="1px solid blue"

}}

InfoA = {{
	width: "100%",
	height: "100%",
	debug: false,
	color: "white",
	addSelectionOptions: false,
	serverURL: "http://chemapps.stolaf.edu/jmol/jsmol/jsmol.php",
	use: use,
  coverImage: "",//"data/1hxw.png",        // initial image instead of applet
  coverScript: "",	// special script for click of cover image (otherwise equal to script)
  deferApplet: false,                  // wait to load applet until click
  deferUncover: false,                 // wait to uncover applet until script completed
  jarPath: "java",
	j2sPath: "j2s",
	jarFile: (useSignedApplet ? "JmolAppletSigned.jar" : "JmolApplet.jar"),
	isSigned: useSignedApplet,
	//disableJ2SLoadMonitor: true,
	//disableInitialConsole: true,
	readyFunction: jmol_isReady,
	//defaultModel: "$dopamine",
	script: "load \\"{}\\" {{2 2 2}};"
}}


}})();



Jmol.setDocument(document);

// notice that we are using no document.write() function here. All DOM-based.
// Jmol.getAppletHtml is working.

$(document).ready(function(){{
  Jmol.setDocument(0);
  jmol1 = Jmol.getApplet("jmol1", InfoA)
  $("#appletplace1").html(Jmol.getAppletHtml(jmol1));
}});

</script>
</head>
<body style="height:100%">
<div id="appletplace1" style="width:80%;height:80%;background-color:red"></div>
</body>
</html>

"""

templates = tornado.template.Loader("web_templates")

# send the index file
class IndexHandler(tornado.web.RequestHandler):
    def get(self, url = '/'):
        self.render('index.html')
    def post(self, url ='/'):
        self.render('index.html')


# handle commands sent from the web browser
class CommandHandler(tornado.web.RequestHandler):
    #both GET and POST requests have the same responses
    def get(self, url = '/'):
        print "get"
        self.handleRequest()

    def post(self, url = '/'):
        print 'post'
        self.handleRequest()

    # handle both GET and POST requests with the same function
    def handleRequest( self ):
        # is op to decide what kind of command is being sent
        op = self.get_argument('op', None)

        #received a "checkup" operation command from the browser:
        if op == "checkup":
            #make a dictionary
            status = {"server": True, "mostRecentSerial": "A LINE"}
            #turn it to JSON and send it to the browser
            self.write(json.dumps(status))

        #operation was not one of the ones that we know how to handle
        else:
            print op
            print self.request
            raise tornado.web.HTTPError(404, "Missing argument 'op' or not recognized")


# handle commands sent from the web browser
class RandomHandler(tornado.web.RequestHandler):
    #both GET and POST requests have the same responses
    def get(self, url = '/'):
        print "get"
        self.handle_request()

    def post(self, url = '/'):
        print 'post'
        self.handle_request()

    # handle both GET and POST requests with the same function
    def handle_request(self):

        # Possible options:
        # replace_only: tuple of sites to replace
        # groups_only: only use specific groups
        # max_different: resrict simultaneous types of groups
        backends = [WebStoreBackend()]
        base_structure = random.choice(initialised_structures)
        random_combination_replace(structure=base_structure, groups=f_groups, backends=backends)

        cif_written = backends[0].cifs.popitem()
        page = templates.load('random.html').generate(cif_name=cif_written[1])
        self.write(page)



# adds event handlers for commands and file requests
application = tornado.web.Application([
    #all commands are sent to http://*:port/com
    #each command is differentiated by the "op" (operation) JSON parameter
    (r"/(com.*)", CommandHandler ),
    (r"/(random.*)", RandomHandler ),
    (r"/", IndexHandler),
    (r"/(index\.html)", tornado.web.StaticFileHandler,{"path": cwd}),
    (r"/(.*\.cif)", tornado.web.StaticFileHandler,{"path": cwd }),
    (r"/(.*\.png)", tornado.web.StaticFileHandler,{"path": cwd }),
    (r"/(.*\.jpg)", tornado.web.StaticFileHandler,{"path": cwd }),
    (r"/(.*\.js)", tornado.web.StaticFileHandler,{"path": cwd }),
    (r"/(.*\.css)", tornado.web.StaticFileHandler,{"path": cwd }),
])


def load_structure(name, data_directory):
    """Load an input structure from the filesystem."""
    # Load an existing pickled structure or generate a new one

    start_directory = os.getcwd()
    os.chdir(data_directory)

    pickle_file = "__%s.fapswitch" % name

    loaded = False

    if path.exists(pickle_file):
        print("Existing structure found: %s; loading..." % pickle_file)
        with open(pickle_file, 'rb') as load_structure:
            input_structure = pickle.load(load_structure)
        # Negative versions ensure that very old caches will be removed
        if not hasattr(input_structure, 'fapswitch_version'):
            input_structure.fapswitch_version = (-1, -1, -1)
        # Need to make sure it is still valid
        if input_structure.fapswitch_version != DOT_FAPSWITCH_VERSION:
            print("Incorrect dot-fapswitch detected, re-initialising")
            loaded = False
        else:
            print("Finished loading")
            loaded = True

    if not loaded:
        print("Initialising a new structure. This may take some time.")
        input_structure = ModifiableStructure(name)
        input_structure.from_file(name, 'cif', None)

        # Ensure that atoms in the structure are properly typed
        input_structure.gen_factional_positions()

        print("Generating topology with Open Babel")
        input_structure.gen_babel_uff_properties()

        # A couple of structure checks
        #input_structure.check_close_contacts()
        #input_structure.bond_length_check()

        # Initialise the sites after bonds are perceived
        input_structure.gen_attachment_sites()
        input_structure.gen_normals()

        # Cache the results
        print("Dumping cache of structure connectivity to %s" % pickle_file)
        print("dot-fapswitch version %r" % (DOT_FAPSWITCH_VERSION, ))
        input_structure.fapswitch_version = DOT_FAPSWITCH_VERSION
        with open(pickle_file, 'wb') as save_structure:
            pickle.dump(input_structure, save_structure)

    os.chdir(start_directory)
    return input_structure





if __name__ == "__main__":

    #
    # Fun fapswitch as a web server
    #
    # Everything lives in the global scope

    # Directory to store pickles

    initial_directory = os.getcwd()

    datastore = os.path.join(os.getcwd(), 'datastore')

    os.chdir(datastore)

    available_structures = glob.glob('*.cif')

    initialised_structures = []

    for available_structure in available_structures:
        sname = available_structure[:-4]
        initialised_structures.append(load_structure(sname, datastore))

    os.chdir(initial_directory)

    # Functional group library is self initialising
    f_groups = FunctionalGroupLibrary()
    print("Groups in library: %s" % str(f_groups.group_list))

    #TODO: set collision tester

    #tell tornado to run checkSerial every 10ms
    #serial_loop = tornado.ioloop.PeriodicCallback(checkSerial, 10)
    #serial_loop.start()

    #start tornado
    application.listen(tornadoPort)
    print("Starting server on port number %i..." % tornadoPort )
    print("Open at http://127.0.0.1:%i/index.html" % tornadoPort )
    tornado.ioloop.IOLoop.instance().start()
