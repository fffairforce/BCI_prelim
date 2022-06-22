import threading
from ..socket import ServiceConsumer

'''
Allows you to enable and disable certain data streams over the network.
'''
class DataClient():
    def __init__(self, server_address):
        self.service = ServiceConsumer('trodes.data.service', server_address=server_address)

    def enable(self, streamname):
        req = {'command': 'enable', 'streamname': streamname}
        self.service.request(req)

    def disable(self, streamname):
        req = {'command': 'disable', 'streamname': streamname}
        self.service.request(req)
