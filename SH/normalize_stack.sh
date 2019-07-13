#!/bin/sh
#
# Abakumov Ivan, Denis Zhurovich
# 29.05.2015 Hamburg University
# e-mail: abakumov_ivan@mail.ru

# Multiply stack on square root of number of traces
# Note: we expect that stack and number of traces have the same size

numtfname='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/test/test15/crsstack.numt.su'

stackfname='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/test/test15/crsstack.stack.su'

resultfname='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/test/test15/result.stack.su'

susort < $numtfname cdp > tmp1.su

segyclean < tmp1.su | sumath op=pow a=0.5 > tmp2.su

susort < $stackfname cdp > tmp3.su

segyclean < tmp3.su > tmp4.su

suprod tmp2.su tmp4.su > $resultfname

rm -f tmp*.su

