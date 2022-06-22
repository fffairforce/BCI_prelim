'''
Represents the XML HardwareConfiguration element 

Includes functions for calculating channels.

Stores no state.
'''
class HardwareChannelMap:
    def __init__(self, channel_map):
        # sort devices based on preference
        self.devices = list(channel_map['HardwareConfiguration']['Device'])
        self.devices.sort(key=lambda x: int(x['_packetOrderPreference']))

    def find_channel(self, channel_name):
        # search channels for matching id
        for device_idx, device in enumerate(self.devices):
            for channel_idx, channel in enumerate(device['Channel']):
                if channel['_id'] == channel_name:
                    return device_idx, channel_idx
        raise KeyError('Could not find channel name in channel map.')
    
    def calculate_digital_index(self, device_idx, channel_idx):
        # offset by the devices coming before
        sync_byte_offset = 1
        device_byte_offset = sum(map(lambda x: int(x['_numBytes']), self.devices[:device_idx]))

        # offset by byte position in device
        channel = self.devices[device_idx]['Channel'][channel_idx]
        channel_byte_offset = int(channel['_startByte'])
        channel_bit = int(channel['_bit'])

        return (sync_byte_offset + device_byte_offset + channel_byte_offset) * 8 + channel_bit
    
    '''
    Assumes 1) fetching interleaved analog data 2) analog channels from previous devices need offset
    '''
    def calculate_analog_index(self, device_idx, channel_idx):
        prev_devices = self.devices[:device_idx]
        prev_devices_channels = [chan for dev in prev_devices for chan in dev['Channel']]
        prev_analog_channels = list(filter(lambda x: self.is_analog(x), prev_devices_channels))

        # calculate the uninterleaved from previous devices
        # each analog channel is 2 bytes
        uninterleaved_offset = len(prev_analog_channels) * 2

        channel = self.devices[device_idx]['Channel'][channel_idx]
        channel_id_byte = uninterleaved_offset + int(channel['_interleavedDataIDByte'])
        channel_id_bit = int(channel['_interleavedDataIDBit'])
        channel_start_byte = uninterleaved_offset + int(channel['_startByte'])

        return channel_id_byte, channel_id_bit, channel_start_byte

    def is_analog(self, chan):
        return chan['_dataType'] == 'analog'