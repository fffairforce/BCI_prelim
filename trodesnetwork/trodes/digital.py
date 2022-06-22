from ..socket import SourceSubscriber
from .hardware_channel_map import HardwareChannelMap

import threading

'''
Use this class to subscribe to analog sources

Requires input of a channel map object. The channel map is just
a JSON-like dictionary of the XML HardwareConfiguration node
in the `.trodesconfig` file.

Requires a server_address to connect to the server.

It can be used like this:

    subscriber = trodes.DigitalClient(
        server_address=self.network_address,
        channel_map=config.channel_map,
        channel_name='ECU_Din8')

'''
class DigitalClient():
    def __init__(self, server_address, channel_map, channel_name):
        self.channel_map = HardwareChannelMap(channel_map)

        device_idx, channel_idx = self.channel_map.find_channel(channel_name)
        index =  self.channel_map.calculate_digital_index(device_idx, channel_idx)

        self.raw_client = RawDigitalClient(server_address=server_address, index=index)

    def receive(self):
        return self.raw_client.receive()

class RawDigitalClient():
    def __init__(self, server_address, index):
        self.index = index
        self.byte_id = index // 8
        self.bit_id = index % 8

        self.subscriber = SourceSubscriber('source.digital', server_address=server_address)

    def receive(self):
        rec = self.subscriber.receive()
        timestamp = rec['localTimestamp']
        data = rec['digitalData'][0]
        bit = (data[self.byte_id] >> self.bit_id) & 1
        return timestamp, bit

'''
Subscriber wraps subscription in a thread and callback

Callback can be used to call a Qt signal
'''
class DigitalSubscriber():
    def __init__(self, server_address, channel_map, channel_name, callback):
        self.thread = DigitalSubscriber.DigitalSubscriberThread(
            server_address, channel_map, channel_name, callback)
        self.thread.start()

    class DigitalSubscriberThread(threading.Thread):
            def __init__(self, server_address, channel_map, channel_name, callback):
                super().__init__(daemon=True)
                self.subscriber = DigitalClient(server_address, channel_map, channel_name)
                self.callback = callback

            def run(self):
                while True:
                    res = self.subscriber.receive()
                    self.callback(res)