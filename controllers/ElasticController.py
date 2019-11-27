#Elastic Controller adjustment by Adeel Khalid SIDDIQUI
import os, sys
import numpy
import time
from utils_py.util import debug, format_bytes, CircularBuffer
from BaseController import BaseController

DEBUG = 1
   
class ElasticController(BaseController):

    def _init_(self):
        super(ElasticController, self)._init_()
        self.queue_target = 10 #40 10 seems to be ok with 20 average is 2.8
        self.k1 = 0.02
        self.k2 = 0.001
        self.t_last = -1
        self.int_error = 0  
        self.prec_state = 0
        self.last_level = 0
        self.bwe_filt = -1
        self.bwe_vec = CircularBuffer(5)
        self.algo = 'Elastic'

    def _repr_(self):
        return '<ElasticController-%d>' %id(self)

    def calcControlAction(self):
        self.bwe_vec.add(self.feedback['bwe'])

    def __harmonic_mean(v):
        x = numpy.array(v)
        return  1.0/(sum(1.0/x)/len(x))
        self.bwe_filt = __harmonic_mean(self.bwe_vec.getBuffer())
        e = self.__getError()
        if e == 0:                             
            self.prec_state = 1
            zero_int_error = 1
        elif e > 0 and self.prec_state == 0: 
            self.prec_state = 2
            zero_int_error = 1
        elif e < 0 and self.prec_state == 2:
            self.prec_state = 0
            zero_int_error = 1
        else: 
            zero_int_error = 0

        max_rate = self.feedback['max_rate']
        min_rate = self.feedback['min_rate']
        d = self.feedback['player_status']
        if self.t_last < 0:
            delta_t = 0
            self.int_error = e
        else:
            delta_t = time.time() - self.t_last
            if zero_int_error == 1:                        
                self.int_error = 0
                self.int_error += delta_t * e

		self.t_last = time.time()
		q = self.feedback['queued_time']
		if q < self._getQueueLowWM() or q > self._getQueueHighWM() or self.feedback['is_check_buffering']:
		    '''The control law is b / ( 1 - k1 e - k2 ei)'''
		    den = 1 - self.k1*e - self.k2*self.int_error
            if self.feedback['is_check_buffering']:
                bwe = self.feedback['bwe']
            else:
                bwe = self.bwe_filt
                u = bwe/den
            if den <= 0 or u >= max_rate:
                u = max_rate + 1000 # 2000 adeel Make sure that the maximum rate can be selected
                debug(DEBUG, '%s calcControlAction: Max rate reached. Anti windup active',self)
                self.__resetIntegralError(delta_t * e)
            elif u <= min_rate:
                u = min_rate 
                self.__resetIntegralError(delta_t * e)
                debug(DEBUG, '%s calcControlAction: Min rate reached. Anti windup active',self)
            else:
                u = self.feedback['cur_rate']

		level_u = self.quantizeRate(u)
		if q > self.__getQueueHighWM() and level_u < self.last_level:
		    u = self.feedback['cur_rate']
		    level_u =  self.last_level
		    debug(DEBUG, '%s Prevented switch down when q > q_H',self)
		self.last_level = level_u
		debug(DEBUG, '%s calcControlAction: u: %s/s q: %.2f e: %.2f int_err: %.2f delta_t: %.2f level_u: %d',  self, 
		    format_bytes(u), q, e, self.int_error, delta_t, level_u)
		return u

    def onPaused(self):
        self.int_error = self.__getError()
        debug(DEBUG, '%s onPaused: Paused. Resetting integral error',self)

    def __resetIntegralError(self, data):
        debug(DEBUG,'%s __resetIntegralError: Resetting integral error',self)
        self.int_error -= data

    def __getQueueLowWM(self):
        return self.queue_target - 1.1*self.feedback['fragment_duration'] #2 adeel

    def __getQueueHighWM(self):
        return self.queue_target + 1.1*self.feedback['fragment_duration'] # 2 adeel

    def __getError(self):
        q = self.feedback['queued_time']
        if q > self.__getQueueHighWM():
            return q - self.__getQueueHighWM()
        elif q < self.__getQueueLowWM():
            return q - self.__getQueueLowWM()
        return 0
