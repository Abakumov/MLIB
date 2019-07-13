#!/bin/sh
#
# Abakumov Ivan
# 12.05.2015 Hamburg University
# e-mail: abakumov_ivan@mail.ru

# Sort output of i-crs code with the cdp parameter


filepath='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/test/test14'

susort < $filepath/icrsstack.coher.su cdp > $filepath/icrsstack.cdp.coher.su
susort < $filepath/icrsstack.numcoher.su cdp > $filepath/icrsstack.cdp.numcoher.su	
susort < $filepath/icrsstack.numt.su cdp > $filepath/icrsstack.cdp.numt.su
susort < $filepath/icrsstack.stack.su cdp > $filepath/icrsstack.cdp.stack.su

#susort < $filepath/crsstack.coher.su cdp > $filepath/crsstack.cdp.coher.su
#susort < $filepath/crsstack.numcoher.su cdp > $filepath/crsstack.cdp.numcoher.su	
#susort < $filepath/crsstack.numt.su cdp > $filepath/crsstack.cdp.numt.su
#susort < $filepath/crsstack.stack.su cdp > $filepath/crsstack.cdp.stack.su


# CMPINI

#filepath='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/test/test9'

#susort < $filepath/cmpini.coher.su     cdp > $filepath/cmpini.cdp.coher.su
#susort < $filepath/cmpini.maxoffset.su cdp > $filepath/cmpini.cdp.maxoffset.su	
#susort < $filepath/cmpini.nip00.su     cdp > $filepath/cmpini.cdp.nip00.su
#susort < $filepath/cmpini.nip10.su     cdp > $filepath/cmpini.cdp.nip10.su
#susort < $filepath/cmpini.nip11.su     cdp > $filepath/cmpini.cdp.nip11.su	
#susort < $filepath/cmpini.numcoher.su  cdp > $filepath/cmpini.cdp.numcoher.su
#susort < $filepath/cmpini.numt.su      cdp > $filepath/cmpini.cdp.numt.su
#susort < $filepath/cmpini.stack.su     cdp > $filepath/cmpini.cdp.stack.su

