from ..socket import ServiceConsumer

'''
'''
class InfoClient:
    def __init__(self, server_address):
        self.service = ServiceConsumer(
            'trodes.info', server_address=server_address)

    def __request(self, item):
        data = { 'request': item }
        return self.service.request(data)

    def request_time(self):
        return self.__request('time')[2]['time']

    def request_timerate(self):
        return self.__request('timerate')[2]['timerate']

    def request_config(self):
        return self.__request('config')