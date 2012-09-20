#!/usr/bin/env python

"""
Fapswitch client for sending functionalisations.

usage: client.py PORT

where PORT is the port that fapswitch is listening on.

"""

# Echo client program
import socket
import sys

# HOST is the machine that fapswitch is running on (i.e. the same machine)
HOST = 'localhost'
# The port depends on which one fapswitch chooses, this just takes the
# second argument to the script.
PORT = int(sys.argv[1])

# Here we create a socket and connect it to the server
# no error checking is done so you wil get an error for
# incorrect ports.
fapswitch = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
fapswitch.connect((HOST, PORT))

# keep going until user types exit
while 1:
    # Get user input from standard input
    line = raw_input("fapswitch >>> ")
    if 'help' in line.lower():
        print("Fapswitch client:")
        print("Freeform strings: dot separated between curly braces")
        print("e.g. {.Me...F..Cl.Cl.Cl..}")
        print("Site replacements: dot separated between square brackets")
        print("e.g. [Me@H7.COOH@H8.F@H12]")
        print("Multiple structures can be input on a single line")
        print("Blank line exits interative mode")
    elif 'exit' in line:
        break
    else:
        # Just send the input data straight to fapswitch
        fapswitch.sendall(line)
        # wait for a reply
        data = fapswitch.recv(1024)
        # dump the reply stright out
        print('fapswitch ::: %s' % data)

fapswitch.close()
