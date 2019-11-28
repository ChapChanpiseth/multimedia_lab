#!/usr/bin/env python
# -*- Mode: Python -*-
# -*- encoding: utf-8 -*-
# Copyright (c) Chap Chanpiseth <chap.chanpiseth@gmail.com>

# This file may be distributed and/or modified under the terms of
# the GNU General Public License version 2 as published by
# the Free Software Foundation.
# This file is distributed without any warranty; without even the implied
# warranty of merchantability or fitness for a particular purpose.
# See "LICENSE" in the source distribution for more information.
import os, sys
from utils_py.util import debug, format_bytes, log
from BaseController import BaseController

from collections import OrderedDict

DEBUG = 1

# This controller is an implementation of the Conventional Controller
# described in Algorithm 1 of the paper:
# Zhi Li, et al, "Probe and Adapt: Rate Adaptation for HTTP Video 
# Streaming At Scale", IEEE JSAC, vol.32, no.4, Apr 2014

class BBA0Controller(BaseController):

    def __init__(self):
        super(BBA0Controller, self).__init__()

        # ---------------------------------------------------
        # BBA-0: Buffer-based
        # ---------------------------------------------------
        self.reservoir = 0.1
        self.cushion = 0.9
        self.buffer_size =60 
        self.initial_buffer = 6
        self.alpha_max= 0.99
        self.delta_t = 10
        
        # Rate map for the bitrates, Rmin to Rmax
        self.rate_map = None
        
        # the algorithm relies on past data, thus as an simple approach we simply save
        # all feedback data in a history list
        self.feedback_hist = []

    def __repr__(self):
        return '<BBA0Controller-%d>' %id(self)

    def isBuffering(self):
        return self.feedback['queued_time'] < self.buffer_size

    def setPlayerFeedback(self, dict_params):
        super(BBA0Controller, self).setPlayerFeedback(dict_params)

        # append the data to our feedback history
        self.feedback_hist.append(self.feedback)

    def calcControlAction(self):
        
        video_rates = self.feedback['rates']
        #curr_bitrate = self.feedback['cur_rate']
        
        # Calculate video rate as next chunk to download
        next_bitrate, Bdelay = self.__adapte_buffer(video_rates, self.rate_map)

        # debug prints
        debug(DEBUG, "%s feedback %s", self, self.feedback)
        # debug(DEBUG, "%s last_download_time %s fragment_duration %s", self, self.feedback['last_download_time'],  self.feedback['fragment_duration'])
        #debug(DEBUG, "%s next_level = %d, Bdelay = %f", self, next_level, Bdelay)
        new_level = self.__levelLessThanRate(next_bitrate)
        
        # The algorithm returns Bdelay which represents the minimum buffer level
        # in seconds of playback when the next download must be started.
        #self.setIdleDuration(0 if Bdelay == 0 else (self.feedback["queued_time"] - Bdelay))
        self.setIdleDuration(0 if Bdelay == 0 else (Bdelay))

        # return next video rate
        return next_bitrate

    def __levelLessThanRate(self, rate):
        vr = self.feedback['rates']
        l = 0
        for i in range(0,len(vr)):
            if rate >= vr[i]:
                l = i
        return l

    def __map_rate(self, video_rates):
        """
            Module to generate the rate map for the bitrates, reservoir, and cushion
        """

        if self.rate_map is None:
            self.rate_map = OrderedDict()
        
        self.rate_map[self.reservoir] = video_rates[0]
        intermediate_levels = video_rates[1:-1]
        marker_length = (self.cushion - self.reservoir)/(len(intermediate_levels)+1)
        current_marker = self.reservoir + marker_length
        for _vid_rate in intermediate_levels:
            self.rate_map[current_marker] = _vid_rate
            current_marker += marker_length
        self.rate_map[self.cushion] = video_rates[-1]

        debug(DEBUG, "%s Map_rate: %s", self, self.rate_map)
        return self.rate_map
    
    # TODO: find a decent explanation
    def time_intersect(self, t1_b, t1_e, t2_b, t2_e):
        if t1_e <= t2_b or t1_b >= t2_e:
            # t2        |----| or |----|
            # t1 |----|                  |----|
            return 0
        elif t1_b < t2_b and t1_e > t2_b:
            # t2    |----|
            # t1 |----|
            return t1_e - t2_b
        elif t1_b < t2_e and t1_e > t2_e:
            # t2 |----|
            # t1    |----|
            return t2_e - t1_b
        else:
            # t2 |----|
            # t1  |--|
            return (t2_e - t2_b) - (t1_b - t2_b) - (t2_e - t1_e)
    
    # calculates the average segment throughput during the time interval [t1, t2]
    def p_tilde(self, t1, t2):
        sum_o = 0
        sum_u = 0

        # loop over every segment
        for feedback in self.feedback_hist[1:]:
            # calculate time difference between segment download start and end
            t_i_b = feedback["start_segment_request"]
            t_i_e = feedback["stop_segment_request"]
            download_duration = t_i_e - t_i_b

            # ignore all segments that have been downloaded after t2
            if t_i_e > t2:
                break

            # get average bit-rate of representation
            p_dach_r_i = feedback["cur_rate"]

            # get segment duration
            tau = feedback["fragment_duration"]

            # calculate segment throughput of segment i
            p_tilde_i = feedback["bwe"]

            sum_o += p_tilde_i * self.time_intersect(t_i_b, t_i_e, t1, t2)
            sum_u += self.time_intersect(t_i_b, t_i_e, t1, t2)

        return sum_o / sum_u
    
    def __adapte_buffer(self, video_rates, rate_map):
        """
        Buffer-based adaptation module
        """
        
        tau = self.feedback["fragment_duration"]
        max_level = self.feedback["max_level"]
        min_level = 0
        current_rate = self.feedback["cur_rate"]
        current_level = self.feedback["level"]
        higher_rate = self.feedback["rates"][current_level if current_level == max_level else (current_level + 1)]
        lower_rate = self.feedback["rates"][current_level if current_level == 0 else (current_level - 1)]
        # debug(DEBUG, "%s Buffer Occupied %s", self, self.feedback["queued_time"])
        buffer_occupied = self.feedback["queued_time"]
        buffer_occupied_pct = buffer_occupied/self.buffer_size
        
        t = self.feedback["stop_segment_request"]
        p_tilde = self.p_tilde(t - self.delta_t, t)
        p_tilde_i = self.feedback["bwe"]
        #debug(DEBUG, "%s ****** throughput calculation %f and current throuput %f and AVg rate %f", self, p_tilde, p_tilde_i,  current_rate)
        
        next_level = current_level
        Bdelay = 0
        lowB = 10
        B_opt = 0.5 * (lowB + self.buffer_size)
        
        # Generate the rate map for the bitrates, Rmin to Rmax
        if not (rate_map):
            rate_map = self.__map_rate(video_rates)
    
        # At the beginning selects the lowest video rate
        if buffer_occupied_pct <= self.reservoir:
            next_vidrate = lower_rate
        elif buffer_occupied_pct >= self.reservoir + self.cushion:
            next_vidrate = higher_rate
        else:
            
            estimated_vidrate =  self.feedback["rates"][current_level]
            
            rate_keys = reversed(rate_map.keys())
            selected_idx = 0.1
            marker = 0.0
            for marker in rate_keys:
                # Find the appropriate rate ratio which quasi matches current buffer rate
                if marker < buffer_occupied_pct:
                    fB = (video_rates[0] * (buffer_occupied - self.buffer_size*self.reservoir))/tau

                    # Estimate video rate
                    if current_level != min_level and current_rate > p_tilde_i:
                        next_level = current_level - 1
                    else:
                       if next_level < max_level: 
                           next_level = current_level + 1
                    
                    debug(DEBUG, "%s ****** Fb %f and Next Vidrate %f", self, fB, self.feedback["rates"][next_level])
                    #if rate_map[marker] <= fB and fB <= self.feedback["rates"][next_level]:
                    if rate_map[marker] <= fB:
                        selected_idx = marker
                        break;

            estimated_vidrate = rate_map[selected_idx]
            debug(DEBUG, "%s Marker %s Estimated_vidrate: %s Throughput %s", self, marker,  estimated_vidrate,  p_tilde_i) 
                
            if (estimated_vidrate >= higher_rate):
                Bdelay = self.buffer_size - tau
                rate_values = reversed(rate_map.values())
                for _vidrate in rate_values:
                    if _vidrate <= higher_rate:
                        next_vidrate = _vidrate
                        debug(DEBUG, "%s Estimated_vidrate > High rate %s", self, next_vidrate)
                        break;
            elif(estimated_vidrate <= lower_rate):
                rate_values = rate_map.values()
                debug(DEBUG, "%s Estimated_vidrate > Lower rate %s estimated_vidrate %s", self, lower_rate, estimated_vidrate)
                for _vidrate in rate_values:
                    if  lower_rate <= _vidrate:
                        next_vidrate = _vidrate
                        debug(DEBUG, "%s Estimated_vidrate > Lower rate %s", self, next_vidrate)
                        break;
            else:
                next_vidrate = current_level
                debug(DEBUG, "%s Estimated_vidrate Default %s", self, next_vidrate)
        
        # Calculate delay before downloading next chunk
        if buffer_occupied>= self.buffer_size:
            T = self.feedback['last_download_time']
            Bdelay =   tau - T
        elif buffer_occupied < self.buffer_size:
            if current_level == max_level or higher_rate  >= self.alpha_max * p_tilde:
                Bdelay = max(buffer_occupied - tau, B_opt)
        else:
            if current_level == max_level or higher_rate >= self.alpha_max * p_tilde:
                    Bdelay = max(buffer_occupied - tau, B_opt)
    
        # return next_vidrate, rate_map, state
        return next_vidrate, Bdelay
