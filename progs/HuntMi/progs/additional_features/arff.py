import sys

header1 = '''@relation dataset

@attribute G+C numeric
@attribute MFEI1 numeric
@attribute MFEI2 numeric
@attribute dG numeric
@attribute dQ numeric
@attribute dF numeric
@attribute ZD numeric
@attribute EAFE numeric
@attribute Diff numeric
@attribute dS numeric
@attribute dS/L numeric
@attribute Avg_Bp_Stem numeric
@attribute |A-U|/L numeric
@attribute |G-C|/L numeric
@attribute |G-U|/L numeric
@attribute (A-U)/stems numeric
@attribute (G-C)/stems numeric
@attribute (G-U)/stems numeric
@attribute MFEI3 numeric
@attribute MFEI4 numeric
@attribute Diversity numeric
@attribute tri1 numeric
@attribute tri2 numeric
@attribute tri3 numeric
@attribute tri4 numeric
@attribute loops numeric
@attribute dustmasker numeric
@attribute translation numeric
@attribute CLASS {0, 1}

@data
'''

header2 = '''@relation dataset

@attribute G+C numeric
@attribute MFEI1 numeric
@attribute MFEI2 numeric
@attribute dG numeric
@attribute dQ numeric
@attribute dF numeric
@attribute ZD numeric
@attribute EAFE numeric
@attribute Diff numeric
@attribute dS numeric
@attribute dS/L numeric
@attribute Avg_Bp_Stem numeric
@attribute |A-U|/L numeric
@attribute |G-C|/L numeric
@attribute |G-U|/L numeric
@attribute (A-U)/stems numeric
@attribute (G-C)/stems numeric
@attribute (G-U)/stems numeric
@attribute MFEI3 numeric
@attribute MFEI4 numeric
@attribute Diversity numeric
@attribute tri1 numeric
@attribute tri2 numeric
@attribute tri3 numeric
@attribute tri4 numeric
@attribute loops numeric
@attribute dustmasker numeric
@attribute translation numeric
@attribute CLASS {0,1}

@data
'''


if len(sys.argv) == 4:

  write = header1

  input_positive = sys.argv[1]
  input_negative = sys.argv[2]
  output = sys.argv[3]

  f_positive = open(input_positive, 'r')
  f_negative = open(input_negative, 'r')
  out = open(output, 'w')

  file_positive = f_positive.readlines()
  file_negative = f_negative.readlines()

  for line in file_positive:
    line = line.replace('\n', '')
    line += ',1' + '\n'
    write += line

  for line in file_negative:
    line = line.replace('\n', '')
    line += ',0' + '\n'
    write += line

if len(sys.argv) == 3:

  write = header2

  input = sys.argv[1]
  output = sys.argv[2]

  f_input = open(input, 'r')
  out = open(output, 'w')

  file_input = f_input.readlines()

  for line in file_input:
    line = line.replace('\n', '')
    line += ',?\n'
    write += line


out.write(write)
  


