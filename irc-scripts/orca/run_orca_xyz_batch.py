#!/usr/bin/env python
from __future__ import print_function
import os
import numpy as np
from os import path, system
from textwrap import dedent
import shutil
import sys # for python version.
import subprocess

#from atom_set import *


if sys.version_info[0] <= 2:
   read_stdin_input = raw_input
else:
   read_stdin_input = input


# this stuff was in atom_set.xyz:
import numpy as np
from numpy import dot, array

ToAng = 0.5291772108  # molpro default.
ToAmu = 1822.88839  # amu in electron masses.
ToCm = 219474.63067
ToKj = 2625.500
ToKcal = 627.5096


def ReadXyzLinesRaw(Lines, Scale):
   # allowed formats: <nAtoms> \n Desc \n <atom-list>
   #              or: <atom-list> (without any headers)
   # in the first case, only the first nAtoms+2 lines are read, in the
   # second case everything which does not look like a xyz line is
   # ignored.
   nAtoms = None
   r = 0,-1
   if ( len(Lines[0].split()) == 1 ):
      nAtoms = int(Lines[0].split()[0])
      Caption = Lines[1]
      r = 2,nAtoms+2
   else:
      Caption = ""
   Atoms = []
   Xyz = []
   for Line in Lines[r[0]:r[1]]:
      ls = Line.split()
      try:
         Atom = ls[0]
         x,y,z = float(ls[1]), float(ls[2]), float(ls[3])
      except:
         continue
      Atom = Atom[0].upper() + Atom[1:].lower()
      # maybe we should allow for (and ignore) group numbers after the
      # elements?
      if Atom not in ElementNames:
         raise Exception("while reading '{}': unrecognized element '{}'.".format(FileName,Atom))
      Atoms.append(Atom)
      Xyz.append((x,y,z))
   Xyz = Scale*array(Xyz).T
   if 0:
      print("*read '{}':\n{}".format(FileName, str(FAtomSet(Xyz, Atoms))))
   return Xyz, Atoms, Caption


def ReadXyzLines(Lines, Scale, Name=None):
   Xyz, Atoms, Caption = ReadXyzLinesRaw(Lines, Scale)
   return FAtomSet(Xyz, Atoms, Name=Name, XyzCaption=Caption)


def ReadXyzFile(FileName,Scale):
   with open(FileName, "r") as File:
      Text = File.read()
   Lines = Text.splitlines()
   return ReadXyzLines(Lines, Scale, Name=FileName)


def ReadMultiXyz(FileName, Reverse=False):
   with open(FileName, "r") as File:
      Lines = File.read().splitlines()
   Energies = []
   Frames = []
   iLine = 0
   while iLine < len(Lines):
      nAtoms = int(Lines[iLine])
      Caption = Lines[iLine+1]
      #for XyzLine in Lines[iLine+2:iLine+2+nAtoms]:
         #xx
      iLineEnd = iLine + 2 + nAtoms
      Frame = ReadXyzLines(Lines[iLine:iLineEnd], Scale=1., Name="{}:{}".format(FileName,len(Frames)))
      Frames.append(Frame)
      Frame.XyzCaption = Caption
      iLine = iLineEnd
   if not Reverse:
      return Frames
   else:
      return Frames[::-1]


def ElementNameDummy():
   ElementNames = "Xx H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn".split()
   ElementNumbers = dict([(o,i) for (i,o) in enumerate(ElementNames)])
   return ElementNames, ElementNumbers
ElementNames, ElementNumbers = ElementNameDummy()


def mdot(*args):
   """chained matrix product: mdot(A,B,C,..) = A*B*C*...
   No attempt is made to optimize the contraction order."""
   r = args[0]
   for a in args[1:]:
      r = dot(r,a)
   return r


def dot2(A,B): return dot(A.flatten(),B.flatten())


class FAtom(object):
   def __init__(self, Element, Position, Index):
      self.Element = Element
      self.Pos = Position
      self.Index = Index
   @property
   def Label(self):
      # return element and center index combined.
      return "{:2}{:3}".format(self.Element,1 + self.Index)
   @property
   def iElement(self):
      return ElementNumbers[self.Element]
   def __str__(self):
      return "{} ({:6.3f},{:6.3f},{:6.3f})".format(self.Label, self.Pos[0], self.Pos[1], self.Pos[2])


class FAtomSet(object):
   def __init__(self, Positions, Elements, Orientations=None, Name=None, XyzCaption=None):
      """Positions: 3 x nAtom matrix. Given in atomic units (ABohr).
      Elements: element name (e.g., H) for each of the positions.
      Orientations: If given, a [3,3,N] array encoding the standard
      orientation of the given atoms (for replicating potentials!). For
      each atom there is a orthogonal 3x3 matrix denoting the ex,ey,ez
      directions."""
      self.Pos = Positions
      assert(self.Pos.shape[0] == 3 and self.Pos.shape[1] == len(Elements))
      self.Elements = Elements
      self.Orientations = Orientations
      self.Name = Name
      self.XyzCaption = XyzCaption  # if from xyz file: stores caption line
   def MakeXyz(self,NumFmt = "%15.8f"):
      Lines = []
      for i in range(len(self.Elements)):
         Lines.append(" %5s {0} {0} {0}".format(NumFmt) % (\
            self.Elements[i], self.Pos[0,i], self.Pos[1,i], self.Pos[2,i]))
      return "\n".join(Lines)
   def nElecNeutral(self):
      """return number of electrons present in the total system if neutral."""
      return sum([ElementNumbers[o] for o in self.Elements])
   def fCoreRepulsion1(self, iAt, jAt):
      if iAt == jAt: return 0. # <- a core doesn't repulse itself.
      ChA, ChB = [ElementNumbers[self.Elements[o]] for o in [iAt, jAt]]
      return ChA * ChB / np.sum((self.Pos[:,iAt] - self.Pos[:,jAt])**2)**.5
   def fCoreRepulsion(self):
      N = len(self.Elements)
      Charges = array([ElementNumbers[o] for o in self.Elements])
      fCoreEnergy = 0
      for i in xrange(N):
         for j in xrange(i):
            fCoreEnergy += self.fCoreRepulsion1(i,j)
            #fCoreEnergy += Charges[i] * Charges[j] / np.sum((self.Pos[:,i] - self.Pos[:,j])**2)**.5
      return fCoreEnergy
   @property
   def XyzEnergy(self):
      # get energy from xyz caption line, if there was one. The format is not standardized.
      # So we just try to return the last thing which looks like a floating point number.
      if self.XyzCaption is None:
         return None
      else:
         L = self.XyzCaption.split()
         for e in L[::-1]:
            try:
               Energy = float(self.XyzCaption[-1])
               return Energy
            except:
               pass
         # if we reached this point there was nothing to return.
         return None


   def __str__(self):         return self.MakeXyz()
   def __len__(self):         return len(self.Elements)
   def __getitem__(self,key): return FAtom(self.Elements[key], self.Pos[:,key], key)
   def __iter__(self):
      for (iAt,(Type,Xyz)) in enumerate(zip(self.Elements, self.Pos.T)):
         #yield (Type,Xyz)
         yield FAtom(Type, Xyz, iAt)




# ---------------------

def GuessEnergyFromXyzCaption(XyzCaption):
   # split up the caption. Search for the first string which looks like a floating point number.
   for s in XyzCaption.split():
      iDot = s.find('.')
      if iDot == -1 or iDot == 0 or iDot == len(s)-1:
         continue
      if not (s[iDot-1].isdigit() and s[iDot+1].isdigit()):
         continue
      try:
         f = float(s)
         return f
      except ValueError as e:
         continue
   return None



class FNebFrame(object):
   def __init__(self, Atoms, Path, iFrame, Options):
      self.Atoms = Atoms
      self.iFrame = iFrame
      self.BaseName1 = "frame{:04d}".format(iFrame)
      #print(self.Atoms.XyzCaption)
      self.Energy = GuessEnergyFromXyzCaption(self.Atoms.XyzCaption)


      if Options.MakeFrameSubPaths:
         # make sub-paths for all individual frames.
         self.BasePath = path.join(Path, self.BaseName1)
      else:
         self.BasePath = Path
      
      self.BaseName = path.join(self.BasePath, self.BaseName1)

   def FileNameForExt(self, Ext):
      return "{}{}".format(self.BaseName, Ext)

   def LoadXml(self):
      return ReadMolproXml(self.FileNameForExt(".xml"), SkipVirtual=True)


class FWfDecl(object):
   def __init__(self, Charge=None, Spin=None, Sym=None, FullDecl=None):
      if FullDecl:
         self.FullDecl = Decl
      else:
         # total charge of system  (+1: Cation, -1: anion)
         self.Charge = Charge
         # spin: 0 -> singlet, 1 -> doublet, 2 -> triplet, etc.
         self.Spin = Spin
         # D2H spatial symmetry irrep. "1" means totally symmetric.
         self.Sym = Sym
         self.FullDecl = None

   def FormatForMolpro(self):
      if self.FullDecl:
         return self.FullDecl
      else:
         L = []
         if self.Charge:
            L.append("charge={}".format(self.Charge))
         if self.Spin:
            L.append("spin={}".format(self.Spin))
         if self.Sym:
            L.append("sym={}".format(self.Sym))
         if L:
            return "wf,{}".format(",".join(L))
         else:
            return ""


def AtomDistance(A, B):
   return np.linalg.norm(A.Pos - B.Pos)

class FImageSeries(object):
   def __init__(self, TargetPath, AtomSets, ComTemplate, WfDecl, Options):
      self.TargetPath = TargetPath
      self.WfDecl = WfDecl
      self.ComTemplate = ComTemplate
      self.Frames = []
      for (iFrame,Frame) in enumerate(AtomSets):
         self.Frames.append(FNebFrame(Frame, TargetPath, iFrame, Options))


   #def RunIboWithOrca(self, OrcaCommand, TmpDir, QsubCmd=None, QsubShTemplate=None, RefGbwFile=None, MoStart=None, iCenterFrame=None):
   def RunIboWithOrca(self, Options):
      Series = self
      if not ((Options.QsubCmd is None) == (Options.QsubShTemplate is None)):
         raise Exception("--qsub_cmd and --qsub_sh_template arguments must either BOTH be given, or neither of them.")
      assert(Options.MoStart in ["none", "ref", "last"])

      iCenterFrame = Options.iCenterFrame
      if iCenterFrame is None:
         iCenterFrame = 0
      elif iCenterFrame.lower() == "ts":
         Energies = [o.Energy for o in Series.Frames]
         if None in Energies:
            raise Exception("sorry, some guessing energy from .xyz failed at least for some frame(s). Can't find maximum energy frame.")
         iCenterFrame = Energies.index(max(Energies))
         print("center_frame=ts -> identified frame #{} with energy {} as center frame. Will trace in both directions from here.".format(iCenterFrame, Energies[iCenterFrame]))
      else:
         iCenterFrame = int(iCenterFrame)

      if iCenterFrame < 0:
         iCenterFrame += len(Series.Frames)
      if iCenterFrame >= len(Series.Frames):
         raise Exception("CenterFrame = {} was specified, but I have only {} frames loaded (note: first frame index is 0!).".format(iCenterFrame, len(Series.Frames)))


      # check if the target path exists. If not, make it.
      if not path.exists(Series.TargetPath):
         print("Target path '{}' does not exist. Create it? (Y/N + Enter)".format(Series.TargetPath))
         Answer = read_stdin_input()
         if Answer.lower() in ["y","yes"]:
            os.makedirs(Series.TargetPath)
            #except OSError,e:
         else:
            raise Exception("Target path '{}' not accessible.", Series.TargetPath)

      assert(Options.TmpDir is None)
      #if Options.TmpDir:
         ## if present, remove an already existing start wf.
         #try:
            #os.remove(path.join(Options.TmpDir, "nebwf.wfu"))
         #except:
            #pass
         ## tell molpro to store all intermediate data in the given temporary directory.
         #MolproCommand = MolproCommand + " -W %(TmpDir)s -I %(TmpDir)s -d %(TmpDir)s" % {"TmpDir": Options.TmpDir}
      CommandList = []
      CommandList.append("#!/bin/bash")
      CommandList.append("set -x #echo on")
      def VarOr0(x):
         if x is None:
            return 0
         else:
            return int(x)
      # FIXME: fix ordering of frames here to optionally start in the middle.
      def RunFrameSubset(OrderedFrames, RefGbwFile):
         StartGbwFile = RefGbwFile
         for Frame in OrderedFrames:
            if Options.MoStart == "none":
               StartGbwFile = None

            if not path.exists(Frame.BasePath):
               os.makedirs(Frame.BasePath)

            ThisGbwFile = Frame.FileNameForExt(".gbw")
            if StartGbwFile is None:
               MoSource = ""
               MoInp = ""
            else:
               MoSource = "moread"
               MoInp = '%moinp "{}"'.format(os.path.relpath(StartGbwFile, Frame.BasePath))

            FmtArgs = {
               "Geometry": Frame.Atoms.MakeXyz(),
               "FileName1": Frame.BaseName1, #without path.
               "FileName": Frame.BaseName,
               "FilePath": Frame.BasePath,
               #"WfDecl": Series.WfDecl.FormatForMolpro(),
               "Charge": VarOr0(Series.WfDecl.Charge),
               "Multiplicity": 1 + VarOr0(Series.WfDecl.Spin),
               "MoSource": MoSource,
               "MoInp": MoInp,
            }
            Input = Series.ComTemplate.format(**FmtArgs)
            InputName = Frame.FileNameForExt(".inp")
            with open(InputName, "w") as File:
               File.write(Input)
            Cmd = "{} {} ".format(Options.OrcaCommand, InputName)
            OrcaPath = path.dirname(Options.OrcaCommand)
            if Options.QsubShTemplate is None:
               # run locally
               CommandList.append("cd '{}'".format(Frame.BasePath))
               CommandList.append(Cmd)
               CommandList.append("{} '{}' -molden".format(os.path.join(OrcaPath, "orca_2mkl"), Frame.BaseName1))
               CommandList.append("")
            else:
               # run on queue
               FmtArgs["OrcaCmd"] = Cmd
               FmtArgs["OrcaPath"] = OrcaPath
               QsubSh = Options.QsubShTemplate.format(**FmtArgs)
               ShName = Frame.FileNameForExt(".sh")
               with open(ShName, "w") as File:
                  File.write(QsubSh)
               os.chmod(ShName, 0o777)
               CommandList.append(Options.QsubCmd.format(Options.QsubCmd, ShName))

            if Options.MoStart == "last":
               StartGbwFile = ThisGbwFile
            #print("!{}".format(Cmd))
            #res = os.system(Cmd)
            #if res != 0:
               #raise Exception("calculation '{}' failed.".format(InputName))
      ForwardFrames = Series.Frames[iCenterFrame:]
      BackwardFrames = Series.Frames[:iCenterFrame][::-1]
      
      RefGbwFileForBackwardFrames = Options.RefGbwFile
      if ForwardFrames:
         RunFrameSubset(ForwardFrames, Options.RefGbwFile)
         RefGbwFileForBackwardFrames = ForwardFrames[0].FileNameForExt(".gbw")
      if BackwardFrames:
         RunFrameSubset(BackwardFrames, RefGbwFileForBackwardFrames)
      
      FullShText = "\n".join(CommandList) + "\n"
      ShFileName = os.path.join(self.TargetPath,'run_all.sh')
      with open(ShFileName, "w") as File:
         File.write(FullShText)
      os.chmod(ShFileName, 0o777)
      print("generated: '{}'".format(ShFileName))


HelpText = """HELP

- Use 'python run_molpro_xyz_batch.py --help' to get a list of supported command
  line arguments.

EXAMPLE INVOKATIONS

- This one makes input files at the target directory (/tmp/neb2), and a file
  /tmp/neb2/run_all.sh containing commands to run all inputs sequentially on
  the local machine:
      python run_molpro_xyz_batch.py ~/Downloads/IRC4.xyz -p /tmp/neb2 --input-template molpro_batch_local.inp --molpro-cmd '/home/cgk/dev/molpro/bin/molpro -n 4'
- this one makes both input files and shell scripts for batch submission:
      python run_molpro_xyz_batch.py ~/Downloads/IRC4.xyz -p /tmp/neb2 --input-template molpro_batch_local.inp --molpro-cmd '/home/cgk/dev/molpro/bin/molpro -n 4' --qsub-cmd 'qsub' --qsub-sh-template molpro_batch_submit.sh

- Commands can be stored in files and loaded. E.g.,
      python run_molpro_xyz_batch.py ~/Downloads/IRC4.xyz @batch.cfg -p /tmp/irc4
  will read the molpro and other commands from file @batch.cfg

COMMAND SUMMARY
"""

def BoolFromText(Desc):
   s = Desc.lower()
   if s in ['y', 'yes', 't', 'true', 'on']:
      return True
   elif s in ['n', 'no', 'f', 'false', 'off']:
      return False
   else:
      raise Exception("Value '{}' not recognized as a truth value. Should be 'yes'/'true' or 'no'/'false'")   


class FInputGenOptions(object):
   def __init__(self, Args):
      self.OrcaCommand = Args.orca_cmd
      if self.OrcaCommand is None:
         try:
            self.OrcaCommand = subprocess.check_output(["which","orca"]).strip()
         except subprocess.CalledProcessError as e:
            raise Exception("\n****\n{}\n****\n^- ERROR: could not locate orca executable via 'which'. Please specify the --orca-cmd command line argument and give the full path to the orca executable".format(str(e), self.OrcaCommand))
         print("Orca executable obtained via 'which orca': '{}'".format(self.OrcaCommand))


      ### -n 10:  10 processor cores
      ### -m 500M:  500 megaword (4000 MB) per processor.
      ##MolproCommand = "$MOLDEV/bin/molpro -n 10"
      ## store temporary wave function and integral data here. Can be set to None. Should be used
      ## if starting guesses are to be transferred between frames.
      self.TmpDir = None
      #if Args.tmp_dir is not None:
         #TmpDir = Args.tmp_dir

      # declare type of wave function (charge, spin)
      self.Wf = FWfDecl(Charge=Args.wf_charge, Spin=Args.wf_spin)

      try:
         with open(Args.input_template, "r") as File:
            self.ComTemplate1 = File.read()
      except os.IOError as e:
         raise Exception("failed to read input template file '{}'".format(Args.input_template))

      ShTemplate1 = None
      if Args.qsub_sh_template is not None:
         try:
            with open(Args.qsub_sh_template, "r") as File:
               ShTemplate1 = File.read()
         except os.IOError as e:
            raise Exception("failed to read queue sh template file '{}'".format(Args.qsub_sh_template))

      if Args.mo_start is not None:
         self.MoStart = Args.mo_start.lower()
      else:
         self.MoStart = "last"
      assert(self.MoStart in ["none", "ref", "last"])

      self.MakeFrameSubPaths = True
      if Args.sub_dirs is not None:
         self.MakeFrameSubPaths = BoolFromText(Args.sub_dirs)
         
      self.QsubCmd = Args.qsub_cmd
      self.QsubShTemplate = ShTemplate1
      self.RefGbwFile = Args.ref_gbw
      self.iCenterFrame = Args.center_frame
      if self.QsubCmd is not None:
         if self.MoStart == "last":
            raise Exception("'--qsub-cmd <...>' is not compatible with '--mo-start last'. Use '--mo-start none' or '--mo-start ref'")
      if self.MoStart == "ref":
         if self.RefGbwFile is None:
            raise Exception("'--mo-start ref' was specified, but no '--ref-gbw <file>' option was given.")
      self.InputNameTemplate = Args.input_name
      if self.InputNameTemplate is None:
         self.InputNameTemplate = "frame{:04d}"
      if '{' not in self.InputNameTemplate or '}' not in self.InputNameTemplate or " " in self.InputNameTemplate.format(1):
         raise Exception("--input-name '{}': The format appears to be broken. A valid example is \"--input-name 'frame{:04d}'\".".format(self.InputNameTemplate))
      pass


def _main():
   import argparse
   ArgParser = argparse.ArgumentParser(fromfile_prefix_chars='@')
   from sys import argv
   import tempfile
   Commands = argv
   if len(argv) == 1:
      print(HelpText)
   if 1:
      # read default arguments from default config file:
      # $HOME/.config/run_orca_xyz_batch.rc
      # if it exists.
      import os
      DefaultConfigFileName = "run_orca_xyz_batch.rc"
      DefaultConfigFileName = os.path.join(os.environ["HOME"], ".config", DefaultConfigFileName)
      if os.path.exists(DefaultConfigFileName):
         Commands = ['@{}'.format(DefaultConfigFileName)] + Commands
   ArgParser.add_argument('-p', '--output-path', required=True, help="path to store input files and output files. Will be created if necessary.")
   #ArgParser.add_argument('-o', '--output-names', required=False, help="template for naming output files, e.g., -o 'ibo%%03i' will make files named ibo000.inp, ibo001.inp, ... By default: Generated from first .xyz input file.")
   ArgParser.add_argument('--input-template', required=True, help="template for ORCA input file. See given example (orca_batch_local.inp). Geometries and wf declarations will be written into this.")
   ArgParser.add_argument('--wf-charge', required=False, help="if given, electronic charge to pass to {Charge} variable in input")
   ArgParser.add_argument('--wf-spin', required=False, help="if given, 1+(this) will be passed as {Multiplicity} variable in input")
   ArgParser.add_argument('--mo-start', required=False, help="one of 'none', 'ref', or 'last' (defaults to 'last' if not given). If 'none', no starting MO commands will be issued. If 'ref', commands for reading MOs (%moinp) from ref-gbw will be issued. If 'last', commands for starting from ref-gbw will be issued for starting frame(s), and for all other frames, gbw's from last frame in respective direction will be taken for '%moinp' ")
   ArgParser.add_argument('--ref-gbw', required=False, help="specify reference .gbw file to use with --mo-start.")
   ArgParser.add_argument('--sub-dirs', required=False, help="one of 'yes' (default), 'no'. If yes, make a sub-directory for the calculation of every frame.")
   ArgParser.add_argument('--center-frame', required=False, help="If given: specify the frame index at which the trace should start. If 'ts', set this as frame of maximum energy. If '0', start at first frame in .xyz and trace forward. If '-1', start at last frame of .xyz file and trace backwards. If a frame index in-between: trace starting from this frame, and proceed into *BOTH* directions.")
   #ArgParser.add_argument('--wf-decl', required=False, help="if given, a wave function declaration to pass to the input scripts as %%(WfDecl)s. E.g., --wf-decl 'wf,charge=1,spin=1'")
   #ArgParser.add_argument('--molpro-cmd', required=True, help="molpro command to execute. '{1}' will be replaced by the base file name (witout extension!). Example: -molpro-cmd '/usr/local/molpro/bin/molpro -n 8 -m 500M --nobackup --no-xml-output {1}.inp'")
   #ArgParser.add_argument('--orca-path', required=True, help="molpro command to execute. Example: -molpro-cmd '/usr/local/molpro/bin/molpro -n 8 -m 500M --nobackup --no-xml-output'")
   ArgParser.add_argument('--orca-cmd', required=False, help="orca executable to use. Example: -orca-cmd '/opt/prog/orca_4_0_1_2_linux_x86-64_shared_openmpi202/orca'. If not given, taken from output of 'which orca'.")
   #ArgParser.add_argument('--tmp-dir', required=False, help="if given, add -W, -I, -d commands with this directory to molpro. Allows accessing files between runs (e.g., restarting orbitals from last frame if running locally)")
   ArgParser.add_argument('--qsub-sh-template', required=False, help="if given, shell script template for scripts to submit to queue submission command. See example molpro_batch_submit.sh")
   ArgParser.add_argument('--qsub-cmd', required=False, help="if given, command to run for queue submission. '{1}' will be replaced shell script name.")
   ArgParser.add_argument('--input-name', required=False, help="describes the name of input files generated. Defaults to 'frame{:04d}'")
   ArgParser.add_argument('input_xyz_files', metavar='N', type=lambda s: str(s), nargs='+', help="one or more .xyz files to read initial geometries from")
   #ArgParser.add_argument('--input-xyz-files', const=[], action='store_const', default=sum)

   # parse command line arguments, possibly joined by the commands from config files (input with
   # leading @).
   Args = ArgParser.parse_args(Commands)
   Options = FInputGenOptions(Args)

   InputXyzFiles = [o for o in Args.input_xyz_files if os.path.splitext(o)[1] != '.py']
   # ^- somehow this contains the script name itself if run with python <...>.py xyz1.xyz xyz2.xyz ...
   print("Input xyz files: {}".format(InputXyzFiles))

   # read input frames from input .xyz files. We will take either single or multi-xyz files for this.
   AtomSets = []
   for InputXyzFile in InputXyzFiles:
      AtomSets += ReadMultiXyz(InputXyzFile)

   # make a image series object which describes how to make input files for the given xyz frames.
   # First argument ("/scr1/reaction_IRCs/claisen") is the output path where the target calculations should go.
   ImageSeries = FImageSeries(Args.output_path, AtomSets, Options.ComTemplate1, Options.Wf, Options)

   #ImageSeries.RunIboWithOrca(OrcaCommand, TmpDir, QsubCmd=Args.qsub_cmd, QsubShTemplate=ShTemplate1, RefGbwFile = Args.ref_gbw, iCenterFrame = Args.center_frame, MoStart = MoStart)
   ImageSeries.RunIboWithOrca(Options)




if __name__ == "__main__":
   _main()
