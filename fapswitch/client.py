#!/usr/bin/env python

"""
Fapswitch client for sending functionalisations.

"""

# Echo client program
import socket
import sys

HOST = 'localhost'
PORT = int(sys.argv[1])
fapswitch = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
fapswitch.connect((HOST, PORT))

while 1:
    line = raw_input("fapswitch >>> ")
    if 'help' in line.lower():
        print("Fapswitch client:")
        print("Random strings: dot separated between curly braces")
        print("e.g. {.Me...F..Cl.Cl.Cl..}")
        print("Site replacements: dot separated between square brackets")
        print("e.g. [Me@H7.COOH@H8.F@H12]")
        print("Multiple structures can be input on a single line")
        print("Blank line exits interative mode")
    elif 'exit' in line:
        break
    else:
        fapswitch.sendall(line)
        data = fapswitch.recv(1024)
        print('fapswitch ::: %s' % data)

fapswitch.close()
