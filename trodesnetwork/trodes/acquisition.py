import threading
from ..socket import ServiceConsumer, SourceSubscriber

'''
Use this class to request changes in the acquisition state.

Requires a server_address to connect to the server.
'''
class AcquisitionClient():
    def __init__(self, server_address):
        self.service = ServiceConsumer('trodes.acquisition.service', server_address=server_address)

    def __request(self, command, time):
        # this can be replaced with objects instead of dicts later
        req = {'command': command, 'timestamp': time}
        self.service.request(req)

    def play(self):
        self.__request('play', 0)

    def pause(self):
        self.__request('pause', 0)

    def stop(self):
        self.__request('stop', 0)

    def seek(self, timestamp):
        self.__request('seek', timestamp)
        pass

'''
Requires callback(timestamp).
'''
class AcquisitionSubscriber():
    def __init__(self, server_address, callback):
        self.server_address = server_address
        self.callback = callback

        self.thread = AcquisitionSubscriber.AcquisitionSubscriberThread(self.server_address, callback)
        self.thread.start()
    
    class AcquisitionSubscriberThread(threading.Thread):
        def __init__(self, server_address, callback):
            super().__init__(daemon=True)
            self.subscriber = SourceSubscriber('trodes.acquisition', server_address=server_address)
            self.callback = callback

        def run(self):
            while True:
                res = self.subscriber.receive()

                if res['command'] == 'seek':
                    timestamp = res['timestamp']
                    self.callback(timestamp)
