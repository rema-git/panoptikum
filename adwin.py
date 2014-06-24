# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 00:07:00 2013

@author: Reinhardt A.W. Maier <rema@zaehlwerk.net>
"""
def GetValue(adwin_parameter, filename='F:\\INBOX\\ADwin\\adwin_code_recent.acm'):
    """
    Parses recent ADwin Code for given parameter value
    """
    import re

    with open(filename, 'rU') as adwin_code:
        regex = re.compile(r'^' + adwin_parameter + '\ =\ (-|)(\d+|\d|\d\.\d|\d+\.\d|\d+\.\d+|\d\.\d+)(|E-\d+|E\d+)$')
        value = None

        for line in adwin_code:
            if len(regex.findall(line)) > 0:
                value = float(line.split()[-1])
                break
        return value

def WriteStandbyCode(filenameOut='F:\\INBOX\\ADwin\\adwin_code_standby.acm', filenameIn='F:\\INBOX\\ADwin\\adwin_code_recent.acm', filenameMask='D:\\Experiment\\Steuerprogramme\\Messprogramme 2013\\LiRb2.0\\standby_mask_mit_agilent.lst'):
    """
    Parses recent ADwin Code and generate a standby file
    """
    with open(filenameMask, 'rb') as mask:
        forbiddenChs = mask.readlines()
        forbiddenChs = [item.rstrip('\r\n') for item in forbiddenChs]

    with open(filenameIn, 'r') as adwin_code:
        with open(filenameOut, 'w') as adwin_code_out:
            for line in adwin_code:
                if any(forbiddenCh in line for forbiddenCh in forbiddenChs):
                    adwin_code_out.write('// STANDBY MODE ---> ' + line)
                else:
                    adwin_code_out.write(line)
