'''
Use this class to query basic information about Trodes.

Requires a network_address to connect to the server.
'''
class TrodesClient():
    def __init__(self, server_address):
        self.server_address = server_address
