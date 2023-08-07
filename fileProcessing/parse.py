def pdbline(l):
    ''' Format the line of text, l, to be able
    to get coords of atoms as float values.
    '''
    # split line of text with empty space as delimiter
    l = list(l.split(' '))
    # vvv remove these strings from the line
    key = lambda x: x not in ('', '\n')
    # filter l by key to remove strings above
    l = list(filter(key,l))
    return l

def str2list(s):
    # for lists of line numbers saved
    # using brackets and commas
    chars = ('[',']','\n',',')
    for z in chars:
        s = s.replace(z,'')
    s = s.split(' ')
    s = filter(lambda x: x!='',s)
    s = tuple(int(x) for x in s)
    return s
