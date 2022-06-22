from ..socket import SinkSubscriber, SourcePublisher, SinkPublisher, SourceSubscriber

import threading

'''
Has an inbox... has a sending thread.

Can only be one on a server.
'''
class EventServer:
    def __init__(self, server_address):
        self.thread = EventServer.EventServerThread(server_address)
        self.thread.start()

    class EventServerThread(threading.Thread):
        def __init__(self, server_address):
            super().__init__(daemon=True)
            self.receiver = SinkSubscriber(
                'trodes.event.inbox',
                server_address=server_address)
            self.publisher = SourcePublisher(
                'trodes.event',
                server_address=server_address)
        
        def run(self):
            while True:
                data = self.receiver.receive()
                self.publisher.publish(data)

class EventSender:
    def __init__(self, server_address, event_name):
        self.event_name = event_name

        self.publisher = SinkPublisher(
            'trodes.event.inbox',
            server_address=server_address)

    def publish(self, localTimestamp, systemTimestamp):
        self.publisher.publish({'name': self.event_name, 'localTimestamp': localTimestamp, 'systemTimestamp': systemTimestamp})

class EventSubscription:
    def __init__(self, server_address, event_name):
        self.event_name = event_name

        self.subscriber = SourceSubscriber(
            'trodes.event',
            server_address=server_address)

    def receive(self):
        res = self.subscriber.receive()
        while res['name'] != self.event_name:
            res = self.subscriber.receive()
        return res

class EventSubscriber:
    def __init__(self, server_address, event_name, callback):
        self.thread = EventSubscriber.EventSubscriberThread(server_address, event_name, callback)
        self.thread.start()

    class EventSubscriberThread(threading.Thread):
        def __init__(self, server_address, event_name, callback):
            super().__init__(daemon=True)
            self.subscriber = EventSubscription(server_address, event_name)
            self.callback = callback
        
        def run(self):
            while True:
                res = self.subscriber.receive()
                self.callback(res)
