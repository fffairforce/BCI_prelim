import threading
from ..socket import ServiceConsumer

'''
Use this class to open and close files from Python.

Requires a server_address to connect to the server.

Future:
- navigate directories
- list directories

'''
class FileClient():
    def __init__(self, server_address):
        self.service = ServiceConsumer('trodes.file.service', server_address=server_address)

    '''
    Opens a file for recording.
    '''
    def open(self, filename):
        # this can be replaced with objects instead of dicts later
        req = {'command': 'open', 'filename': filename}
        self.service.request(req)

    '''
    Starts recording on the open file
    '''
    def start(self):
        # this can be replaced with objects instead of dicts later
        req = {'command': 'start', 'filename': ''}
        self.service.request(req)

    '''
    Pauses recording on the open file
    '''
    def pause(self):
        # this can be replaced with objects instead of dicts later
        req = {'command': 'pause', 'filename': ''}
        self.service.request(req)

    '''
    Closes a file for recording. A file is assumed to already be open.
    '''
    def close(self):
        # this can be replaced with objects instead of dicts later
        req = {'command': 'close', 'filename': ''}
        self.service.request(req)

class FileStatusSubscriber():
    pass