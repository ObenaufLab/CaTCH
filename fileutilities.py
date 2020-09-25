#!/usr/bin/env python3

## version CaTCH 0.7.1


import os, sys, string, re, subprocess, random, argparse
import pandas as pd
from builtins import list
from collections import Counter


#####   F U N C T I O N S   #####


# http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
def natural_sorted(l):
    """Sort list of numbers/strings in human-friendly order.

    Args:
        l(list): A list of strings.
    Returns:
        list
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def expand_fpaths(flist):
        """Fully expand and absolute-ify the paths of listed files.

        Does not verify path validity. All paths are expanded.

        Args:
            flist[str]: A list/FilesList of files.
        Returns:
            [str]: List of expanded paths.
        """
        return [os.path.abspath(os.path.expanduser(str(f))) for f in flist]


# Helper function.
def autonumerate(things):
    """Detect duplicate entries in a string list and suffix them.

    Suffixes are in _N format where N a natural number >=2. Existing suffixes
    in that format will also be detected and incremented.

    Args:
        things[str]: A list of strings.
    Returns:
        [str]: A corrected list of strings.
    """
    c = Counter(things);
    # Because I use decrement, reversing the list ensures first instance gets smallest number.
    things.reverse()
    for i, t in enumerate(things):
        n = c[t]
        if n > 1:  # The first occurrence is not suffixed.
            newname = t +'_' + str(n)
            while newname in things:  # Check for already present suffixes
                n += 1
                newname = t +'_' + str(n)
            things[i] = newname
            c[t] -= 1
    things.reverse()
    return things


def make_names(items, parameters):
    """Automatically create file names based on parameters.

    If automatic names happen to turn out identical with one another, unique
    numbers are appended to differentiate them. Check documentation for
    autonumerate().

    Args:
        items[str]: A list of strings/filenames/paths to use as the basis for
                    the output names.
        parameters(str,str,str): The first element is the output directory,
                    the second is a common prefix to add to the names,
                    the third is a common suffix to add to the names.
                    Like so: <out[0]>/<out[1]>item<out[2] .
                    If any of the 3 values in None, no outnames will be made.
                    Use current directory and empty strings as necessary.
    Returns:
        [str]: A list of file paths.
    """
    outfiles = []
    if None not in parameters:
        for i in items:
            outfiles.append(os.path.join(os.path.abspath(os.path.expanduser(parameters[0])),
                                         parameters[1] + i + parameters[2]) )
        autonumerate(outfiles)
    return outfiles


def do_foreach(flist, comm, progress=True, out=(None,None,None), log=False):
    """Execute an arbitrary command for each of the listed files.

    Enables executing a shell command over a range of items, by inserting the
    item values into the command as directed by place-holder substrings.
    Although the above is how it is meant to be used, the values in the
    FilesList could be overridden to be any arbitrary string, in which case,
    only {val} will have the desired effect. The other placeholder values are
    computed with the assumption of the values being files, so may not be
    sensible when the items are not files.

    This is the only function with comments or progress attributes, because
    invoked commands print their own output directly, so any informative messages
    controlled by this library will need to be inserted in real time.

    Args:
        flist[]: A FilesList.
        comm[str]: The components of an arbitrary command, with place-holders:
                    {abs} : absolute path of file.
                    {dir} : absolute path of the file's directory.
                    {val} : the actual value specified as target
                    {bas} : the basename of the file, without the last extension.
                    {cor} : the basename of the file, without any extensions.
                    {ali} : the alias for the file, if iterating through a FilesList.
                    Placeholders can be nested, to allow nested calls of fileutilities:
                    i.e. {{abs}}. A layer of nesting is peeled off each time the function is called,
                    until the placeholders are allowed to be evaluated.
        progress(bool): Show start and completion of iterations on STDERR.
                    (Default True)
        out(str,str,str): The first element is the output directory, the second
                    is a common prefix to add to the names, the third is a
                    common suffix to add to the names. Check documentation for
                    make_names().
        log(bool): Log to /commands.log each individual call.
    """
    outstream= sys.stdout
    # Create output files. [] if out contains None.
    outfiles = make_names(flist, out)
    for i, (myfile, myalias) in flist.enum():
        # Substitute place-holders.
        command = []
        for c in comm:
            # Evaluate placeholders, if they are not nested.
            (mypath, mybase) = os.path.split(str(myfile))
            c = re.sub(r"(?<!\{){abs}(?!\})", str(myfile), c)                   # absolute full path to file
            c = re.sub(r"(?<!\{){dir}(?!\})", mypath, c)                        # absolute path of directory
            c = re.sub(r"(?<!\{){val}(?!\})", mybase, c)                        # filename
            c = re.sub(r"(?<!\{){bas}(?!\})", os.path.splitext(mybase)[0], c)   # filename minus last extension
            c = re.sub(r"(?<!\{){cor}(?!\})", mybase.split('.')[0], c)          # filename minus all extensions
            c = re.sub(r"(?<!\{){ali}(?!\})", str(myalias), c)                  # custom or automatic alias
            # Peel off a layer of nesting for the remaining placeholders and flags.
            c = c.replace('{{abs}}', '{abs}')
            c = c.replace('{{dir}}', '{dir}')
            c = c.replace('{{val}}', '{val}')
            c = c.replace('{{bas}}', '{bas}')
            c = c.replace('{{cor}}', '{cor}')
            c = c.replace('{{ali}}', '{ali}')
            c = c.replace(',-', '-')
            # This argument is ready to go now.
            command.append(c)
        # Redirect output.
        if outfiles:
            outstream = open(outfiles[i], 'w')
        # Verbose stuff.
        see = " ".join(command)
        # Do the thing.
        subprocess.call(" ".join(command), stdout=outstream, shell=True)
        # More verbose stuff.
        if outfiles:
            outstream.close()



#####   C L A S S E S   #####



class FilesList(list):
    """A container for a list of files.

    An extension of the built-in list, specifically for files, providing a
    means to import multiple filenames either from text lists or from
    directories. The purpose is to facilitate batch operations and sensible
    output of their results.

    FilesList is generally backwards compatible with the built-in list and it
    should be possible for them to be used interchangeably in most cases. A
    plain list can be cast as a FilesList, when necessary, allowing appointment
    of default alias values. A FilesList should always work as a plain list
    without additional actions (except for its string representation). When a
    FilesList is accessed as a plain list, only the full paths will be
    accessed. Certain equivalent methods are supplied for

    Most basic operations inherited from list are supported. Appending has been
    overridden to keep paths and aliases in sync. Sorting, deleting and
    inserting of items are currently not supported and will break the
    correspondence between paths and aliases.

    Attributes defined here:
        aliases = [] : Practical aliases for the full file-paths.
    """
    def __init__(self, files=None, aliases=None, fromtuples=None, verbatim=True):
        """Construct an instance of the FilesList.

        A FilesList can be created:
        - empty
        - from a list of files (with default aliases automatically assigned)
        - from a list of files and a list of aliases (in the same order)
        - from a list of (file, alias) tuples.

        Args:
            verbatim(bool): Whether to skip path pre-preocessing. (Default True)
            files[str]: A list of files. (Default None)
            aliases[str]: A list of aliases. (Default None)
            fromtuples[(str,str)]: A list of tuples (file, alias). (Default
                        None) If this is specified together with flist and/or
                        aliases, the data in fromtuples is used only.
        """
        # If data is a list of (file, alias) tuples, unpair tuples into lists.
        if fromtuples is not None:
            data = [list(t) for t in zip(*fromtuples)]
            # Any data passed to flist and aliases at method call is discarded.
            files = data[0]
            aliases = data[1]
        # Having aliases without the paths is rather useless.
        if aliases:
            if not files:
                raise ValueError("No files supplied for the aliases.")
        else:
            # Make None into empty.
            aliases = []
        # Assign default aliases to be same as files. Expand file paths.
        if files is not None:
            if not verbatim:
                files = expand_fpaths(files)
            if not aliases:
                for f in files:
                    aliases.append(self.autoalias(f))
        else:
            # If still empty, it was an empty call to the constructor.
            files = []
        # Create the basic list.
        super().__init__(files)
        # Add a plain list attribute for the aliases with default values.
        self.aliases = autonumerate(aliases)

    def __str__(self):
        """Represent as string.

        Overrides built-in list's representation.

        Returns:
            str
        """
        tmp = []
        for f, (myfile, myalias) in self.enum():
            tmp.append("\t".join([str(f), myfile, myalias]))
        tmp.append("")
        return "\n".join(tmp)

    def to_file(self, outfile=None, mode ='a'):
        """Save list as a text file that can be read back in.

        Args:
            outfile(str): Output file to write into. If omitted, it only
                        returns the content as a print-ready string.
                        (Default None)
            mode(str): Append ('a') or overwrite ('w'). (Default 'a')
        Returns:
            str: A print-ready multi-line string. This is returned even when an
                        output file is specified and written into.
        """
        result = ""
        for f, (myfile, myalias) in self.enum():
            result += myfile + "\t" + myalias + "\n"
        if outfile is not None:
            with open(outfile, mode) as out:
                out.write(result)
        return result

    def enum(self):
        """Enumerate as (index, (filepath, filealias)).

        Returns:
            enumerator"""
        return enumerate(zip(self, self.aliases))

    def get(self, loc):
        """Access path and alias at specified location as tuple.

        Args:
            loc[int]: Index of item to get.
        Returns:
            (str,str): Tuple (path, alias).
        """
        return (self[loc], self.aliases[loc])

    def append(self, myfile, myalias=None, verbatim=True):
        """Appends value to both the paths list and the aliases list.

        This method overrides the built-in append() of list. It is backwards
        compatible by automatically guessing an alias.
        This reduces the risk of the paths and aliases going out of sync due
        to items being manually added without updating the aliases.
        It is still possible to break the sync by manually adding items to the
        aliases.

        Args:
            myfile(str): File (path will be expanded).
            myalias(str): Alias for the file (Default None).
            verbatim(bool): Do not pre-process path for the target value. (Default True)
        """
        if myfile is not None:
            if not verbatim:
                myfile = expand_fpaths([myfile])[0]
        super().append(myfile)
        if not myalias:
            myalias = self.autoalias(myfile)
        self.aliases.append(myalias)
        self.aliases = autonumerate(self.aliases)

    def populate_from_files(self, myfiles, colSep="\t", verbatim=True, alias_verbatim=True):
        """Parse the list of files from one or multiple text files.

        Read in multiple lists of files from text and append them to the
        FilesList. All paths are automatically expanded and converted to
        absolute paths. Because of this, each file may also be given a more
        convenient alias. If no alias is given, the filename as supplied is
        used as the alias instead. The paths are not checked for existence.

        Existing contents of the object are kept and the new contents are
        appended.

        Input file format (no spaces allowed inside names):

            #comment
            path1/file1     alias1-prefix    alias1-suffix1     alias1-suffix2
            path1/file2     alias1-prefix    alias1-suffix3
            path2/file3     alias3
            path3/file4

        Args:
            file[str]: A list of text files each containing a list of files.
            colSep(str): Column separator. (Default "\\t")
            verbatim(bool): Do not pre-process paths for the target values. (Default True)
        Returns:
            FilesList: Returns self, to facilitate instantiation shortcuts.
        """
        # Read new list.
        paths = []
        for myfile in myfiles:
            with open(myfile, 'rU') as input:
                for line in input:
                    if line == "\n":
                        # Skip empty lines.
                        continue
                    elif line[0] == '#':
                        # Skip comments.
                        continue
                    else:
                        fields = line.rstrip().split(colSep)
                        paths.append(fields[0])
                        # Store the alias for the file.
                        if len(fields) > 1:
                            self.aliases.append("_".join(fields[1:]))
                        # If an alias was not specified, re-use the filepath given.
                        else:
                            self.aliases.append(self.autoalias(fields[0]))
        # Expand to absolute paths and add to main self list.
        if not verbatim:
            paths = expand_fpaths(paths)
        self.extend(paths)
        if not alias_verbatim:
            self.aliases = autonumerate(self.aliases)
        return self

    def populate_from_directories(self, dirpaths, patterns=None, verbatim=True, alias_verbatim=True):
        """Get files based on naming patterns from within a list of directories.

        Useful for selecting among files that follow a naming convention. The
        convention is represented by a list of regex strings, at least one of
        which has to match.
        File paths will be expanded automatically. The filenames will be used
        as the aliases.

        Existing contents of the object are kept and the new contents are
        appended.

        Args:
            dirpaths[str]: A list/FilesList of paths to directories from where
                        to get files.
            patterns[str]: A list of regex strings. Only files matching at least
                        one of these will be returned. The patterns will be
                        matched anywhere in the filenames.
            verbatim(bool): Do not pre-process paths for the target values. (Default True)
        Returns:
            FilesList: Returns self, to facilitate instantiation shortcuts.
        """
        rx = []
        if patterns:
            rx = [re.compile(p) for p in patterns]
        for d in dirpaths:
            try:
                os.chdir(d)
                for f in os.listdir(d):
                    if f in ["","\n",".",".."]:
                        continue
                    if not patterns:
                        # No filter.
                        self.append(os.path.join(d, f), self.autoalias(f), verbatim=verbatim)
                    else:
                        for p in rx:
                            if p.search(f):
                                self.append(os.path.join(d, f), self.autoalias(f), verbatim=verbatim)
                                break
            finally:
                pass
        self.aliases = autonumerate(self.aliases)
        return self.sorted()

    # Helper function.
    @staticmethod
    def autoalias(pathname):
        """Strip a path to the base filename."""
        if pathname is None:
            return None
        else:
            return os.path.splitext(os.path.basename(str(pathname)))[0]

    def sorted(self):
        """Sorted copy.

        Returns:
            FilesList
        """
        d = dict()
        for i, (myfile, myalias) in self.enum():
            d[myfile] = myalias
        sk = natural_sorted(list(d.keys()))
        newFL = FilesList()
        for k in sk:
            newFL.append(k, d[k])
        return newFL




#####   M A I N   #####


def main(args):
    """Trimmed version of my fileutilities.py library/script.

    Provided to ensure independence from the continuous development and status of my full script.
    Removed all functionality not relevant to its use in CaTCH.
    """

    # Organize arguments and usage help:
    parser = argparse.ArgumentParser(description="Provide INPUTTYPE and TARGETs *before* providing any of the other parameters. This is due to many parameters accepting an indefinite number of values. Only one task at a time.")

    # Input/Output.
    parser.add_argument('INPUTTYPE', type=str, choices=['T','P'],
                                help=" Specify the type of the TARGETs: 'T' = The actual input filess. 'L' = Text file(s) listing the input files. 'P' = Get list of input files from STDIN pipe. 'D' = Input data directly from STDIN pipe. ('D' is compatible with only some of the functions)")
    parser.add_argument('TARGET', type=str, nargs='*',
                                help=" The targets, space- or comma-separated. Usually files. Look into the specific task details below for special uses. Do not specify with INPUTTYPE 'P' or 'D'.")
    parser.add_argument('--dir', type=str, nargs='*',
                                help=" List the contents of the target paths. Full absolute file paths are returned. Each file is also given an alias. Supplying an optional list of regex patterns enables filtering of the result.")
    parser.add_argument('--loop', type=str, nargs='+',
                                help=" Repeat the specified shell command for each target value. Available PLACEHOLDERS to insert the targets into the commands: {abs} full path, {dir} path of directory portion, {val} target value such as filename, {bas} basename (filename minus outermost extension), {cor} filename core (all extensions removed), {ali} file alias. Flags intended for the nested command should be preceded by a ',' like this: ',-v'. Recursive calls to fileutilities.py are possible by nesting the placeholders and escapes: i.e. {{abs}}, ,,-v. One layer is peeled off with each call to fileutilities loop. The placeholders will take the values of the targets of the respectively nested call.")
    params = parser.parse_args(args)

    # INPUT ###################################################################

    targets = []
    for t in params.TARGET:
        v = t.split(",")
        if len(v) == 1:
            targets.append(t)
        else:
            targets.extend(v)

    flist = None
    if params.INPUTTYPE == 'P':
        # Read files list from STDIN
        flist = FilesList()
        for line in sys.stdin:
            fields = line.rstrip("\n").split("\t")
            if fields[0] != "":
                try:
                    flist.append(fields[0], fields[1])
                except IndexError:
                    flist.append(fields[0])
    elif params.INPUTTYPE == 'T':
        # Create the FilesList by supplying a direct list of files.
        flist = FilesList(targets, verbatim=False)
    else:
        sys.exit("Unknown INPUTTYPE.")

    outdir, outpref, outsuff = None, None, None

    # TASKS ###################################################################

    # Filter DIRECTORY contents. ----------------------------------------------
    if params.dir is not None:
        result = FilesList().populate_from_directories(flist, params.dir)
        sys.stdout.write(result.to_file())

    # LOOP a command. -------------------------------------------------
    elif params.loop:
        command = []
        for c in params.loop:
            command.append(c.lstrip("+"))
        do_foreach(flist, command, out=(outdir, outpref, outsuff),
                   progress=(False), log=False)



#####    E X E C U T I O N   #####


# Call main only if the module was executed directly.
if __name__ == "__main__":
    main(sys.argv[1:])

    sys.exit(0)

#EOF
