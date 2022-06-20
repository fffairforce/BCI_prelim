from trodesnetwork import trodes
from trodesnetwork import socket as trodeSocket
from trodesnetwork.socket import SinkSubscriber
from trodesnetwork.socket import SourcePublisher
from trodesnetwork.socket import SinkPublisher
from trodesnetwork.socket import SourceSubscriber

import socket as s

TRODE_DATA_TYPE = 'source.lfp'
TRODE_ADDRESS = 'tcp://127.0.0.1'
TRODE_PORT = ':49152'
TRODE_SERVER_ADDRESS = TRODE_ADDRESS+TRODE_PORT

UDP_IP = "127.0.0.1"
UDP_PORT = 5005
MESSAGE = b" "

sock = s.socket(s.AF_INET, # Internet
                     s.SOCK_DGRAM) # UDP

lfpSub = trodeSocket.SourceSubscriber(TRODE_DATA_TYPE,server_address=TRODE_SERVER_ADDRESS)
print("Begin Main Loop")
while(True):

    data = lfpSub.receive();
    MESSAGE=str(data).encode('ascii')
    sock.sendto(MESSAGE, (UDP_IP, UDP_PORT))
