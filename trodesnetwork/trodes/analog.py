import threading
from ..socket import SourceSubscriber
from .hardware_channel_map import HardwareChannelMap
import struct

'''
Use this class to subscribe to analog sources

`channel_map` is just a dictionary of the channels XML config.
'''
class AnalogSubscription():
    def __init__(self, server_address, channel_map, channel_name):
        self.channel_map = HardwareChannelMap(channel_map)

        device_idx, channel_idx = self.channel_map.find_channel(channel_name)
        id_byte, id_bit, start_byte = self.channel_map.calculate_analog_index(device_idx, channel_idx)

        self.raw_subscriber = RawAnalogSubscription(
            server_address=server_address,
            id_byte=id_byte,
            id_bit=id_bit,
            start_byte=start_byte)

    def receive(self):
        return self.raw_subscriber.receive()

'''
index is the 0-based index of the analog channels 
'''
class RawAnalogSubscription():
    def __init__(self, server_address, id_byte, id_bit, start_byte):
        self.id_byte = id_byte
        self.id_bit = id_bit
        self.start_byte = start_byte

        self.subscriber = SourceSubscriber('source.analog', server_address=server_address)

    def receive(self):
        rec = self.subscriber.receive()
        # merge segments together
        data = bytes().join(rec['analogData'])

        if data[self.id_byte] & (1 << self.id_bit):
            value, = struct.unpack('>h', data[self.start_byte:self.start_byte+2])
            timestamp = rec['localTimestamp']
            return timestamp, value
        else:
            # make this iterative, not recursive
            return self.receive()

'''
Subscriber wraps subscription in a thread and callback

Callback can be used to call a Qt signal
'''
class AnalogSubscriber():
    def __init__(self, server_address, channel_map, channel_name, callback):
        self.thread = AnalogSubscriber.AnalogSubscriberThread(
            server_address, channel_map, channel_name, callback)
        self.thread.start()

    class AnalogSubscriberThread(threading.Thread):
            def __init__(self, server_address, channel_map, channel_name, callback):
                super().__init__(daemon=True)
                self.subscriber = AnalogSubscription(server_address, channel_map, channel_name)
                self.callback = callback

            def run(self):
                while True:
                    res = self.subscriber.receive()
                    self.callback(res)
