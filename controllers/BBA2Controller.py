#!/usr/bin/env python
# -*- Mode: Python -*-
# -*- encoding: utf-8 -*-

import os, sys
from utils_py.util import debug, format_bytes, log
from BaseController import BaseController
from collections import OrderedDict

DEBUG = 1

# This controller is an implementation of the Conventional Controller
# described in Algorithm 1 of the paper:
# Zhi Li, et al, "Probe and Adapt: Rate Adaptation for HTTP Video 
# Streaming At Scale", IEEE JSAC, vol.32, no.4, Apr 2014

class BBA2Controller(BaseController):

    def __init__(self):
        super(BBA2Controller, self).__init__()

        # ---------------------------------------------------
        # BBA-2: Buffer-based
        # ---------------------------------------------------
        self.reservoir = 0.1
        self.cushion = 0.9
        self.buffer_size = 60 #40
        self.initial_buffer = 2
        self.alpha_max= 0.99
        self.delta_t = 7
        self.initial_factor = 0.875
        self.alpha_4 = 0.75
        self.delta_t = 10
        
        # Rate map for the bitrates, Rmin to Rmax
        self.rate_map = None
        
        # the algorithm relies on past data, thus as an simple approach we simply save
        # all feedback data in a history list
        self.feedback_hist = []
        
        # Start up state
        self.state = "INITIAL"
        
        # shows whether we are in fast start mode
        self.runningFastStart = True
        
        # algorithm configuration parameters
        self.conf = {
            "B_min": 10,
            "B_low": 10,
            "B_high": 60,
            "delta_beta": 1,
            "delta_t": 10,
            "alpha_1": 0.75,
            "alpha_2": 0.33,
            "alpha_3": 0.5,
            "alpha_4": 0.75,
            "alpha_5": 0.9
        }

    def __repr__(self):
        return '<BBA2Controller-%d>' %id(self)

    def isBuffering(self):
        return self.feedback['queued_time'] < self.buffer_size

    def setPlayerFeedback(self, dict_params):
        super(BBA2Controller, self).setPlayerFeedback(dict_params)

        # append the data to our feedback history
        self.feedback_hist.append(self.feedback)

    def calcControlAction(self):
        
        video_rates = self.feedback['rates']
        #curr_bitrate = self.feedback['cur_rate']
        
        # Calculate video rate as next chunk to download
        next_bitrate, Bdelay = self.__adapte_buffer(video_rates, self.rate_map)

        # The algorithm returns Bdelay which represents the minimum buffer level
        # in seconds of playback when the next download must be started.
        if self.feedback["queued_time"] < Bdelay:
            Bdelay = Bdelay - self.feedback["queued_time"]
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

            # ignore all segments that have been downloaded after t2
            if t_i_e > t2:
                break

            # calculate segment throughput of segment i
            p_tilde_i = feedback["bwe"]

            sum_o += p_tilde_i * self.time_intersect(t_i_b, t_i_e, t1, t2)
            sum_u += self.time_intersect(t_i_b, t_i_e, t1, t2)

        return sum_o / sum_u
    
    def __next_vidrate(self,  buffer_occupied_pct,  rate_map,  current_rate,  current_level,  tau,  min_level,  max_level,   B_opt):
        """ Next Video Rate
        """
        
        Bdelay = 0
        next_level = self.feedback["level"]
        p_tilde_i = self.feedback["bwe"]
        higher_rate = self.feedback["rates"][current_level if current_level == max_level else (current_level + 1)]
        
        estimated_vidrate = current_rate
        
        # At the beginning selects the lowest video rate
        if buffer_occupied_pct >= self.reservoir + self.cushion:
            next_vidrate = higher_rate
        else:
            estimated_vidrate =  self.feedback["rates"][current_level]
            rate_keys = reversed(rate_map.keys())
            selected_idx = 0.1
            marker = 0.0
            for marker in rate_keys:
                # Find the appropriate rate ratio which quasi matches current buffer rate
                # Estimate video rate
                if current_level != min_level and current_rate > p_tilde_i:
                    next_level = current_level - 1
                else:
                   if next_level < max_level: 
                       next_level = current_level + 1
                       #next_level = max_level
            estimated_vidrate = self.feedback["rates"][next_level]

        next_vidrate = estimated_vidrate
        # Return next video rate
        Bdelay = 0
        return next_vidrate,  Bdelay

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
        
        next_level = current_level
        Bdelay = 0
        lowB = 10
        B_opt = 0.3 * (lowB + self.buffer_size)
        
        next_vidrate = self.feedback["rates"][current_level]
        
        # Generate the rate map for the bitrates, Rmin to Rmax
        if not (rate_map):
            rate_map = self.__map_rate(video_rates)
            next_vidrate = lower_rate
        elif self.state == "INITIAL":
            # if the B increases by more than 0.875V s. Since B = V - ChunkSize/c[k],
            # B > 0:875V also means that the chunk is downloaded eight times faster than it is played
            next_vidrate = current_rate
            # delta-B = V - ChunkSize/c[k]
            delta_B = tau - current_rate/p_tilde_i            
            # Select the higher bitrate as long as delta B > 0.875 * V
            if delta_B > self.initial_factor * tau:
                # Next video rate
                next_vidrate = self.feedback["rates"][next_level + 1]
    
            # if the current buffer occupancy is less that NETFLIX_INITIAL_BUFFER, then do NOY use rate map
            if buffer_occupied_pct > self.reservoir:
                for _idx in range(max_level+1):
                    _rate = self.feedback["rates"][_idx ]
                    if (_rate > p_tilde_i):
                        next_level = _idx - 1
                        break
    
                # Next video rate
                est_next_vidrate = self.feedback["rates"][next_level]
                
                # Once the rate map returns a higher value exit the 'INITIAL' stage
                if est_next_vidrate >= next_vidrate:
                    next_vidrate = est_next_vidrate
                    T = self.feedback['last_download_time']
                    Bdelay = Bdelay - self.feedback["queued_time"]
                    self.state = "RUNNING"
    
                if current_level == max_level or next_vidrate  >= self.alpha_4 * p_tilde:
                    Bdelay = max(buffer_occupied - tau, B_opt)
                    Bdelay = Bdelay - self.feedback["queued_time"]
        else:
            # At the beginning selects the lowest video rate
            if buffer_occupied_pct <= self.reservoir:
                # self.state = "INITIAL"
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
                        
                        #if rate_map[marker] <= fB and fB <= self.feedback["rates"][next_level]:
                        if rate_map[marker] <= fB:
                            selected_idx = marker
                            break

                estimated_vidrate = rate_map[selected_idx]
                    
                if (estimated_vidrate >= higher_rate):
                    Bdelay = self.buffer_size - tau
                    rate_values = reversed(rate_map.values())
                    for _vidrate in rate_values:
                        if _vidrate <= higher_rate:
                            next_vidrate = _vidrate
                            break
                elif(estimated_vidrate <= lower_rate):
                    rate_values = rate_map.values()
                    for _vidrate in rate_values:
                        if  lower_rate <= _vidrate:
                            next_vidrate = _vidrate
                            break
                else:
                    next_vidrate = current_level
        
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
