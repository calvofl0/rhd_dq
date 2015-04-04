#!/usr/bin/env
# -*- coding: utf-8 -*-

from threading import Thread
from SimpleXMLRPCServer import SimpleXMLRPCServer
import xmlrpclib
import cPickle as pickle
import zlib

class Server(Thread):
	def __init__(self, host='localhost', port=3000):
		Thread.__init__(self)
		self._host = host
		self._port = port
		self._server = SimpleXMLRPCServer((self._host,self._port))
		self._running = False
		self._dict = None
		self.start()
	def set(self,dictionary):
		self._dict = dictionary
	def run(self):
		self._server.register_function(lambda:True,'ping')
		def get_var(var_str):
			d = self._dict
			exec 'q = d.'+var_str in locals()
			var=q
			if hasattr(var, 'value'): var = var.value
			cstr = zlib.compress(pickle.dumps(var, pickle.HIGHEST_PROTOCOL),9)
			return xmlrpclib.Binary(cstr)
		self._server.register_function(get_var,'get')
		self._running = True
		while self._running:
			self._server.handle_request()
		self._server.server_close()
	def stop(self):
		self._running = False
		proxy = xmlrpclib.ServerProxy('http://'+self._host+':'+str(self._port)+'/')
		proxy.ping()

class Client(object):
	def __init__(self, host='localhost', port=3000):
		self._host = host
		self._port = port
		self._proxy = xmlrpclib.ServerProxy('http://'+self._host+':'+str(self._port)+'/')
	def get(self, var_str):
		stream = self._proxy.get(var_str)
		return pickle.loads(zlib.decompress(stream.data))
