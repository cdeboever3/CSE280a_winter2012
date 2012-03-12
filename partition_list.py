def partition_list(lst, parts):
    """
    Partitions the list <lst> into <parts> sub_lists of almost equal
    length.  The first <parts> - 1 sublists will be of equal length,
    and the last sublist will take the reminder.  For example, if
    lst=[1,2,3,4,5] and parts = 3, then [[1,2],[3,4],[5]] is returned.
    """

    # The length of the first parts -1 sublists
    section = int(math.ceil(len(lst) / float(parts)))
    
    partitions = []
    for i in range(parts):
        partitions.append(lst[section*i:min(len(lst), section*(i+1))])
    return partitions

