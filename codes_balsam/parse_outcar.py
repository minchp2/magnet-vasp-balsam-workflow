from collections import defaultdict

NIONS = 0

def parse_equals(line, offset, delim='='):
    return line[offset:].split(delim)[1].split()[0]

def find_parse_equals(line, keyword, delim='='):
    return line[line.find(keyword):].split(delim)[1].split()[0]

def parse_table(lines, line_no, values_offset=2, header_offset=2, nrows=None, strip_text="# of ion", ncols=None):
    global NIONS
    if nrows is None:
        nrows = NIONS
    header = lines[line_no+header_offset].strip().strip(strip_text).split()
    if ncols is None:
        ncols = len(header)
    start = line_no+header_offset+values_offset
    end = start+nrows
    values = [list(map(float, l.split()[-ncols:])) for l in lines[start:end]]
    return header, values

def parse_magnetization(line,offset,line_no,lines):
    idx = line.find("(") + 1
    component = line[idx]
    header, values = parse_table(lines, line_no)
    return (component, header, values)

def parse_totalcharge(line,offset,line_no,lines):
    if line.strip() != 'total charge': return
    return parse_table(lines, line_no)

def parse_stress(line,offset,line_no,lines):
    return parse_table(lines, line_no, header_offset=1, 
            values_offset=12, nrows=1, strip_text="Direction"
    )

def parse_latvecs(line,offset,line_no,lines):
    header, values =  parse_table(
        lines,line_no, header_offset=0,
        values_offset=1,nrows=3,ncols=6,
    )
    return values

def parse_nions(line,offset,line_no,lines):
    global NIONS
    nions = int(parse_equals(line,offset))
    NIONS = nions
    return nions

def parse_pos_force(line,offset,line_no,lines):
    header, values = parse_table(
        lines, line_no,
        header_offset=0,
        values_offset=2,
        ncols=6,
    )
    return values


rules = {
    'NBANDS': lambda line,offset,line_no,lines: int(parse_equals(line,offset)),
    'NKDIM': lambda line,offset,line_no,lines: int(parse_equals(line,offset)),
    'NIONS': parse_nions,
    'FREE ENERGIE OF THE ION-ELECTRON SYSTEM': 
        lambda line,offset,line_no,lines: float(find_parse_equals(lines[line_no+2],'TOTEN')),
    'LOOP+' : lambda line,offset,line_no,lines: float(line.split()[-1]),
    'magnetization (' : parse_magnetization,
    'total charge' : parse_totalcharge,
    'FORCE on cell' : parse_stress,
    'volume of cell' : lambda line,offset,line_no,lines: float(parse_equals(line,offset,delim=':')),
    'direct lattice vectors' : parse_latvecs,
    'TOTAL-FORCE' : parse_pos_force,
}

labels = {
    'FREE ENERGIE OF THE ION-ELECTRON SYSTEM' : "energy",
    'magnetization (' : 'magnetization',
    'FORCE on cell' : 'stress',
    'direct lattice vectors' : 'lattice_vecs',
    'TOTAL-FORCE' : 'position_force',
}


def parse_outcar():
    result = defaultdict(list)
    with open('OUTCAR') as fp:
        lines = fp.readlines()

    for line_no, line in enumerate(lines):
        for keyword, extractor in rules.items():
            offset = line.find(keyword)
            if offset >= 0:
                try:
                    dat = extractor(line, offset, line_no, lines)
                except Exception as e:
                    print(f"[line #{line_no}] Failed to parse {keyword} starting from:\n{line}")
                else:
                    if dat is not None: result[keyword].append(dat)
    result = dict(result)
    for key in list(result.keys()):
        if len(result[key]) == 1: result[key] = result[key][0]
        if key in labels:
            newkey = labels[key]
            result[newkey] = result[key]
            del result[key]
    return result

if __name__ == "__main__":
    from pprint import pprint
    res = parse_outcar()
    for key,val in res.items():
        print(key)
        if isinstance(val, list):
            pprint(val[-1])
        else:
            pprint(val)
