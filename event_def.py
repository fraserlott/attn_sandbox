#!/usr/local/sci/bin/python2.7
'''
@author: Andrew Ciavarella
Module to contain functions relevant to the definition and handling of events.
Contains Event class.
'''

__author__ = "Andrew Ciavarella"

class Event:
    def __init__(self, framing, p0, p1):
        # Framing information
        var = framing[0]
        period = framing[1]
        reg = framing[2]
        yc = framing[3]
        # Varible short name, e.g. 'tas', 'pr'
        self.var = var
        # Averaging period, e.g. 'month', 'day'
        self.period = period
        # Region info
        self.reg = reg
        # Threshold, if threshold, could be None
        self.yc = yc

        # Probabilities as lists of three values: [min, best, max]
        self.p0 = p0
        self.p1 = p1
        try:#if len(p0)==3:
            far_best = 1.0 - (p0[1]/p1[1])
            far_min = 1.0 - (p0[2]/p1[0])
            far_max = 1.0 - (p0[0]/p1[2])
            self.far = [far_min, far_best, far_max]
        except IndexError:#elif len(p0)==1:
            self.far = 1.0 - (p0/p1)
        except TypeError:
            self.far = 1.0 - (p0/p1)
        
    def define(self):
        print('--- Event: %s %s %s' % (self.var, self.period, self.reg))
        # State values to 2 significant figures
        print('    Threshold, yc = %.2g' % self.yc)
        print('    p1 = %.2g [%.2g, %.2g]' % (self.p1[1], self.p1[0], self.p1[2]))
        print('    p0 = %.2g [%.2g, %.2g]' % (self.p0[1], self.p0[0], self.p0[2]))
        print('    FAR = %.2g [%.2g, %.2g]' % (self.far[1], self.far[0], self.far[2]))

def get_test_events(n):
    '''
    Return list of event objects with synthetic test data. 
    '''
    import numpy as np
    # Probabilities ultimately to be passed to function
    p1s = [(np.random.random())**2 for i in range(n)]
    p1s_min = [np.random.uniform(0.0,p1s[i]) for i in range(n)]
    p1s_max = [np.random.uniform(p1s[i],1.0) for i in range(n)]
    p0s = [np.random.uniform(0,p1s[i]) for i in range(n)]
    p0s_min = [np.random.uniform(0.0,p0s[i]) for i in range(n)]
    p0s_max = [np.random.uniform(p0s[i],1.0) for i in range(n)]   
    # Framing data
    var = 'tas'
    period = 'month'
    reg = ['Reg %s' % (i+1) for i in range(n)]
    yc = [np.random.normal(2.0,1.0) for i in range(n)]
    events = []
    for i in range(n):
        p0_event = [p0s_min[i], p0s[i], p0s_max[i]]
        p1_event = [p1s_min[i], p1s[i], p1s_max[i]]
        f = [var, period, reg[i], yc[i]]
        x = Event(f, p0_event, p1_event)
        #x.define()
        events.append(x)

    return events

########################################################
####################### MAIN ###########################
########################################################

def main():
    '''
    @author: Andrew Ciavarella
    Module to contain functions relevant to the definition and handling of events.
    Contains Event class.
    '''

if __name__ == "__main__":
    main()


