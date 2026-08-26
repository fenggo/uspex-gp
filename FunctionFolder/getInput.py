#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
USPEX 9.4.4 release
2015 Oganov's Lab. All rights reserved.
"""
from __future__ import print_function
import sys
import os

from optparse import OptionParser

# -------------------------------- Options --------------------------------------
parser = OptionParser()
parser.add_option("-f", "--input_txt", dest="input_txt",
                  help="Specify INPUT.txt file", metavar="INPUT")
parser.add_option("-b", "--begin_keyword", dest="begin_keyword",
                  help="Specify begin keyword", metavar="BEGIN")
parser.add_option("-e", "--end_keyword", dest="end_keyword",
                  help="Specify end keyword", metavar="END")
parser.add_option("-c", "--col_number", dest="col_number",
                  help="Specify column number", metavar="COL")

(options, args) = parser.parse_args()


if not options.input_txt:
   parser.error('INPUT.txt file not specified.')
elif not options.begin_keyword:
   parser.error('Keyword/keyblock not specified.')

if os.path.abspath(options.input_txt):
   input_txt = options.input_txt
else:
    parser.error('Input file ' + options.input_txt + ' not found.')

begin_keyword = options.begin_keyword
end_keyword = None
col_number = None

if options.end_keyword:
    end_keyword = options.end_keyword
else:
    # Column number must be specified:
    if not options.col_number:
        parser.error('Column number not specified.')
    else:
        col_number = options.col_number

# print('-  ',col_number)
# parser.destroy()


# Read INPUT.txt file
f = open(input_txt, 'rb')
input_content = f.readlines()
f.close()

print('<CALLRESULT>')

# -------------------------------------------------------------------------------
# Process the case of keyblock:
begin_num = None
end_num = None
if end_keyword is not None:
    for i in range(len(input_content)):
        begin_find = str(input_content[i]).lower().find(begin_keyword.lower())
        if begin_find >= 0 and begin_find < 10 and not begin_num:
            begin_num = i + 1
        line = str(input_content[i].lower())
        if line.find(end_keyword.lower()) >= 0 and not end_num:
           end_num = i - 1
           break

    if begin_num is not None and end_num is not None:
       for i in range(begin_num, end_num + 1):
           key = str(input_content[i].strip())
           # print(key,key[0],key[1])
           if key[0] == 'b':
              key = key[2:]
              key = key[:-1]
           print(key)

# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------

str_num = None
if end_keyword is None:
   for i in range(len(input_content)):
       line = str(input_content[i].lower())
       if line.find(begin_keyword.lower()) >= 0:
          str_num = i
          break

if str_num is not None:
   line=str(input_content[str_num])
   key = line.split(':')[0].split()[int(col_number) - 1]
   if key[0] == 'b':
      key = key[2:]
   print(key)
# -------------------------------------------------------------------------------

sys.exit(0)
