from setuptools import setup

setup(
    name="trodesnetwork",
    version="0.0.11",
    description="A library to connect to Trodes over a network",
    packages=["trodesnetwork", "trodesnetwork.trodes"],
    install_requires=[
        'pyzmq >=18.0.0,<20.0.0',
        'msgpack'
    ]
)

