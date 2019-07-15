% SEG2_READ - Read in a file in SEG-2 format
% Version 0.9
% Usage:
%
%  [FDStrings,nTraces,TDStrings,TraceData]=seg2_read(fname)
%     read data and header from file 'fname'.
%  
%  Outputs:
%     FDStrings - cell array of strings in the File Description Header
%     nTraces   - number of traces read
%     TDStrings - cell array of cell arrays of the Trace Description strings
%                 (e.g., TDStrings{2}{3} is the third string in the second trace)
%     TraceData - cell array of column vectors containing the data
%                 (e.g., TraceData{3} is the third trace).
%
%  Cell arrays are used in the output to allow different trace lengths.
%  Note: because of widespread non-compliance to the SEG-2 standard in
%  field data no attempt is made to extract meaning from the header strings.
%  This is left to a higher level program, probably on a case by case basis.
%
%  (C) James Wookey, University of Bristol, 2008.
%  Reference: SEG-2 standard: http://www.seg.org/publications/tech-stand/seg_2.pdf

%  This software is distributed under the term of the BSD free software license.
%
%  Copyright:
%     (c) 2003-2008, James Wookey, University of Bristol
%
%  All rights reserved.
%
%   * Redistribution and use in source and binary forms, with or without
%     modification, are permitted provided that the following conditions are
%     met:
%        
%   * Redistributions of source code must retain the above copyright notice,
%     this list of conditions and the following disclaimer.
%        
%   * Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.
%     
%   * Neither the name of the copyright holder nor the names of its
%     contributors may be used to endorse or promote products derived from
%     this software without specific prior written permission.
%
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
%   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
%   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
%   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [FDStrings,nTraces,TDStrings,TraceData]=seg2_read(fname) ;
   if nargout~=4, error('Wrong number of output arguments supplied.'), end ;
   fid=fopen(fname,'r') ;

%  read identifier block
   endianFormat='l' ;

   byteCount = 0 ;
   
   id=fread(fid,1,'int16',0,endianFormat) ;
   
   if id==14933
%  ** format is little endian (do nothing)
   elseif id==21818
%  ** format is big endian
      endianFormat='b' ;
   else
%  ** identifier is invalid
      error('Bad identifier at SEG-2 file start, can''t determine format');
   end
   
   revno=fread(fid,1,'int16',0,endianFormat) ;
   SizeTPSubBlock=fread(fid,1,'int16',0,endianFormat);
   if rem(SizeTPSubBlock,4)~=0, error('Bad size for TPSubBlock'); ,end;
   
   nTP = SizeTPSubBlock/4 ;
   nTraces=fread(fid,1,'int16',0,endianFormat) ;
   
   i=fread(fid,1,'int8',0,endianFormat) ; 
   code=fread(fid,[1 2],'int8',0,endianFormat) ;
   strTermCode = code(1:i) ; 
   strTermChar = char(code(1:i)) ;
   i=fread(fid,1,'int8',0,endianFormat) ; 
   code=fread(fid,[1 2],'int8',0,endianFormat) ;
   lineTermCode = code(1:i) ;
   lineTermChar = char(code(1:i)) ; 

%% skip the unused part of the header   
   dump = fread(fid,18,'int8',0,endianFormat) ;

   byteCount = byteCount+32;

   TP_tmp=fread(fid,[1 nTP],'uint32',0,'l') ;
   TP=TP_tmp(1:nTraces) ;

   byteCount = byteCount+nTP*4;
   
   
   nchar=TP(1)-32-nTP*4; % first offset minus the header

   % ~order 806
   code=fread(fid,[1 nchar],'int8',0,endianFormat) ; 

   byteCount = byteCount+nchar ;

%  parse the string structure in here. 
   [FDStrings]=getstrings(code,strTermCode) ;
   
   TDStrings = cell(1,nTraces);
   TDnsamp = cell(1,nTraces) ;
   TraceData = cell(1,nTraces) ;
   
   for ii=1:nTraces

   trace_header_id = fread(fid,1,'uint16',0,endianFormat) ;
   if trace_header_id~=17442
      error('Non-standards compliant SEG-2 format, cannot continue') ;
   end
   
   trace_header_size = fread(fid,1,'uint16',0,endianFormat)  ;
   
   trace_size = fread(fid,1,'uint32',0,endianFormat) ;
   nsamp = fread(fid,1,'uint32',0,endianFormat) ;
   TDnsamp{ii} = nsamp ;
   
   format_code = fread(fid,1,'uint8',0,endianFormat) ;
   
   if format_code~=4
      error('Sorry, only 32 floating point data is currently supported')
   end
   
   dump = fread(fid,19,'uint8',0,endianFormat) ;   


   code=fread(fid,[1 trace_header_size-32],'int8',0,endianFormat);  
%  parse the string structure in here. 
   [TDStrings{ii}]=getstrings(code,strTermCode); 

   byteCount = byteCount+trace_header_size;
   
   trace=fread(fid,[1 nsamp],'float32',0,endianFormat) ;

   TraceData{ii}=trace ;
   
   byteCount = byteCount + nsamp*4;
   
   end
   
   
   fclose(fid);
return

function [STR]=getstrings(code,strTermCode)   
%  first use the string term char in order to split
   strEndsR=find(code==strTermCode); 
   strStartsR=[1 strEndsR(1:(length(strEndsR)-1))+1];
   strEndsR = strEndsR - 1;
   strLengthsR = strEndsR - strStartsR;
%  dump any zero length strings
   [strLengths ind]=find(strLengthsR>0);
   strStarts = strStartsR(ind); 
   strEnds = strEndsR(ind);
   
   nstr = length(strStarts);
   STR = cell(1,nstr) ;
   for istr=1:nstr
      STR{istr} = char(code(strStarts(istr):strEnds(istr)));
   end   


