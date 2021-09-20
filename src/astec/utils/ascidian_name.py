

def get_daughter_names(name):
    #
    # build daughter names from parent name
    #
    # anterior or posterior character 'a' or 'b'
    # stage (round of division)
    # '.'
    # p value (cell index)
    # left or right character '_' or '*'
    #

    # fake returns
    if name == 'background':
        return ['background', 'background']
    if name == 'other-half':
        return ['other-half', 'other-half']

    abvalue = name.split('.')[0][0]
    stage = name.split('.')[0][1:]
    p = name.split('.')[1][0:4]
    lrvalue = name.split('.')[1][4]
    #
    # build daughter names
    #
    daughters = [abvalue + str(int(stage) + 1) + "." + '{:0{width}d}'.format(2 * int(p) - 1, width=4) + lrvalue,
                 abvalue + str(int(stage) + 1) + "." + '{:0{width}d}'.format(2 * int(p), width=4) + lrvalue]
    # print("name = " + str(name) + " -> daughter names = " + str(daughters))
    return daughters


def get_mother_name(name):
    #
    # build daughter names from parent name
    #
    # anterior or posterior character 'a' or 'b'
    # stage (round of division)
    # '.'
    # p value (cell index)
    # left or right character '_' or '*'
    #

    # fake returns
    if name == 'background':
        return 'background'
    if name == 'other-half':
        return 'other-half'

    abvalue = name.split('.')[0][0]
    stage = name.split('.')[0][1:]
    p = name.split('.')[1][0:4]
    lrvalue = name.split('.')[1][4]
    #
    # build parent names
    #
    if int(stage) == 0:
        return None
    parent = abvalue + str(int(stage) - 1) + "."
    if int(p) % 2 == 1:
        parent += '{:0{width}d}'.format((int(p) + 1) // 2, width=4)
    else:
        parent += '{:0{width}d}'.format(int(p) // 2, width=4)
    parent += lrvalue
    # print("name = " + str(name) + " -> parent name = " + str(parent))
    return parent


def get_sister_name(name):
    sister_names = get_daughter_names(get_mother_name(name))
    sister_names.remove(name)
    return sister_names[0]


def get_symmetric_name(name):
    symname = name[:-1]
    if name[-1] == '*':
        symname += '_'
    elif name[-1] == '_':
        symname += '*'
    else:
        return None
    return symname
