from pathlib import Path
import datetime
import csv

def _read_csv(filepath, delimiter=',', comments='#', strip=True):
    """Return a dictionary 

    Args:
        filepath (string or Path): csv filepath.
        delimiter (string, optional): Delimiter string. Default is ','
        comments (string, optional): string that signals comment lines.
            Default is '#'.
        strip (bool, optional): If True, trailling spaces in each item are 
            removed. In practice, this allows for column aligned csv files.
            Default is True.

    Returns:
        dict
    """
    with open(Path(filepath)) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=delimiter)
        values = []
        for row in csv_reader:
            if len(row) > 0:  # remove empty lines
                if row[0].startswith(comments) == False: # remove comments
                    if strip:
                        values.append([_.strip() for _ in row])
                    else:
                        values.append(row)

    return {_[0]: _[1:] for _ in values}

def _str2datetime(string):
    """convert I21 date string pattern to date --> '2022-07-20T21:08:36.921'

    Args:
        string (str): string with I21 date string

    Returns:
        datetime.datetime object
    """
    #########
    # split #
    #########
    date, time = string.split('T')
    
    ########
    # date #
    ########
    year, month,  day = (int(_) for _ in date.split('-'))

    ########
    # time #
    ########
    hour, minute, seconds = time.split(':')

    hour    = int(hour)
    minute  = int(minute)
    seconds = int(round(float(seconds)))

    if seconds >= 60:
        minute  = minute + 1
        seconds = 0
    if minute >= 60:
        hour   = hour + 1
        minute = 0

    ############
    # datetime #
    ############
    return datetime.datetime(year=year, month=month, day=day, hour=hour, minute=minute, second=seconds)