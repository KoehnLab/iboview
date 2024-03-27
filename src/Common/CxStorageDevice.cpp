/* Copyright (c) 2015  Gerald Knizia
 * 
 * This file is part of the IboView program (see: http://www.iboview.org)
 * 
 * IboView is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * IboView is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IboView (LICENSE). If not, see http://www.gnu.org/licenses/
 * 
 * Please see IboView documentation in README.txt for:
 * -- A list of included external software and their licenses. The included
 *    external software's copyright is not touched by this agreement.
 * -- Notes on re-distribution and contributions to/further development of
 *    the IboView software
 */

#include <stdexcept>
#include <stdio.h> // for tmpfile and c-style file handling functions
#include <algorithm> // for std::min
#include <memory.h>
#include <string.h>
#include <stdlib.h>
#include <time.h> // for srand for temp file names... yes, I know how bad both of those are (seeding with time and the cstdlib rand in particular)

#include <fcntl.h>
#include <errno.h>

#include "CxStorageDevice.h"
#include "CxMemoryStack.h"
#include "CxAlgebra.h"
#include "format.h"

namespace ct {

FStorageBlock FStorageDevice::AllocNewBlock(FStreamOffset SizeInBytes)
{
    FRecord r = AllocNewRecord(SizeInBytes);
    return FStorageBlock(*this, r);
}

FStorageDevice::~FStorageDevice()
{
}



FRecord FStorageDevicePosixFs::AllocNewRecord(FStreamOffset SizeInBytes)
{
    FILE *NewFile = tmpfile();
    if (NewFile == 0)
        throw std::runtime_error("FStorageDevicePosixFs: Failed to open a temporary file via tmpfile() function.");
    return AllocNewRecord(SizeInBytes, NewFile);
}

FRecord FStorageDevicePosixFs::AllocNewRecord(FStreamOffset SizeInBytes, FILE* FileHandle)
{
    uint
        NewId = m_FileIds.size();
    if (FileHandle == 0)
        throw std::runtime_error("FStorageDevicePosixFs::AllocNewRecord: Input file handle was NULL.");
    if (m_ControlledFiles.find(FileHandle) != m_ControlledFiles.end())
        throw std::runtime_error("FStorageDevicePosixFs::AllocNewRecord: An input file handle was supplied for multiple different records. This is not allowed.");
    m_FileIds[NewId] = FileHandle;
    m_ControlledFiles.insert(FileHandle);
    return FRecord(NewId, 0, 0, SizeInBytes);
}

FRecord FStorageDevicePosixFs::AllocNewRecord_TempFileInDir(FStreamOffset SizeInBytes, std::string const &Directory, std::string const &FileNameTemplate)
{
    // find a file name which was not yet used.
    bool srand_called = false;
    if (!srand_called) {
        srand(time(NULL));
    }
    FILE
        *NewFile = 0;
    std::string
        FileNameAbsolute,
        FileNameBase = FileNameTemplate,
        RandomPart(10, '-');
    char const
        Letters[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
    size_t const
        nLetters = (sizeof(Letters)/sizeof(Letters[0]) - 1); // -1: it counts the 0-terminator....
    for (size_t iTries = 0; iTries != 0xfffff; ++ iTries) {
        // make a semi-random string to put into our new file name.
        for (size_t i = 0; i != RandomPart.size(); ++ i)
            RandomPart[i] = Letters[size_t(rand()) % nLetters];
        // assemble the total file name, including directory (if given).
        FileNameBase = fmt::format(FileNameTemplate, RandomPart);
        if (Directory.empty()) {
            FileNameAbsolute = FileNameBase;
        } else {
            if (Directory[Directory.size()-1] == '/')
                FileNameAbsolute = Directory + FileNameBase;
            else
                FileNameAbsolute = fmt::format("{}/{}", Directory, FileNameBase);
        }
        // try to open the file for writing, provided it does not already exist.
        // Unfortunately, in C++98 there seems to be no race-condition free way of doing this.
        // C2011 has mode 'w+bx' which fails if the file already exists, but this is not currently
        // part of the C++ standard. POSIX also has some modes which should work for this:
        // https://stackoverflow.com/questions/230062/whats-the-best-way-to-check-if-a-file-exists-in-c-cross-platform
        // ...but this would make the program unaccessible in Windows.
        // So at the moment: Just try and ...hope for the best.
        NewFile = fopen(FileNameAbsolute.c_str(), "w+bx");
        if (NewFile != 0) {
            // at least in linux this will keep the file alive as long as
            // the handle is open...  and since the "w+bx" (With the 'x')
            // would most likely anyway have already crashed the program by
            // now if this wasn't linux, maybe that is okay?
            remove(FileNameAbsolute.c_str());
            break;
        }
    }
    if (NewFile == 0)
        throw std::runtime_error(fmt::format("FStorageDevicePosixFs::AllocNewRecord_TempFileInDir: Failed to create temporary file (dir='{}', filename='{}').", Directory, FileNameBase));
    return AllocNewRecord(SizeInBytes, NewFile);
}



void FStorageDevicePosixFs::Write(FRecord const &r, void const *pData, FStreamOffset nLength, FStreamOffset Offset) const
{
    FILE *File = GetHandle(r);
    fseek(File, Offset + r.BaseOffset, SEEK_SET);
    size_t nWritten = fwrite(pData, 1, nLength, File);
    if (nWritten != nLength || ferror(File)) throw std::runtime_error(fmt::format("FStorageDevicePosixFs:: failed to write {} bytes of data at offset {}.", nLength, Offset));
}

void FStorageDevicePosixFs::Read(FRecord const &r, void *pData, FStreamOffset nLength, FStreamOffset Offset) const
{
    FILE *File = GetHandle(r);
    fseek(File, Offset + r.BaseOffset, SEEK_SET);
    size_t nRead = fread(pData, 1, nLength, File);
    if (nRead != nLength || ferror(File)) throw std::runtime_error(fmt::format("FStorageDevicePosixFs:: failed to read {} bytes of data at offset {}.", nLength, Offset));
}

void FStorageDevicePosixFs::Reserve(FRecord const &/*r*/, FStreamOffset /*nLength*/) const
{
}

void FStorageDevicePosixFs::Delete(FRecord const &r)
{
    FILE *File = GetHandle(r);
    fclose(File);
    m_ControlledFiles.erase(File);
    m_FileIds.erase(r.iFile);
}

FILE *FStorageDevicePosixFs::GetHandle(FRecord const &r) const {
    FFileIdMap::const_iterator
        it = m_FileIds.find(r.iFile);
    if (it == m_FileIds.end())
        throw std::runtime_error("FStorageDevicePosixFs: Attempted to access a non-exitent file id.");
    return it->second;
}


// FILE *FStorageDevicePosixFs::GetHandle(FRecord const &r) {
//     return const_cast<FILE*>( const_cast<FStorageDevicePosixFs const*>(this)->GetHandle(r) );
// };


FStorageDevicePosixFs::FStorageDevicePosixFs()
{
}


FStorageDevicePosixFs::~FStorageDevicePosixFs()
{
    FFileIdMap::reverse_iterator
        itFile;
    for (itFile = m_FileIds.rbegin(); itFile != m_FileIds.rend(); ++ itFile)
        fclose(itFile->second);
}


void FStorageDevice::FRecord::SetLength(FStreamOffset NewLength) const
{
    assert(NewLength != UnknownSize);
    assert(EndOffset == UnknownSize || EndOffset - BaseOffset == NewLength);
    const_cast<FStorageDevice::FRecord*>(this)->EndOffset = BaseOffset + NewLength;
}








FStorageDeviceMemoryBuf::FBuffer &FStorageDeviceMemoryBuf::GetBuf(FRecord const &r) const {
    assert(r.iRecord < m_Buffers.size());
    return *(const_cast<FBuffer*>(&m_Buffers[r.iRecord]));
}

FRecord FStorageDeviceMemoryBuf::AllocNewRecord(FStreamOffset SizeInBytes)
{
    m_Buffers.push_back(FBuffer());
    m_Buffers.back().p = 0;
    m_Buffers.back().Length = 0;
    FRecord
        r = FRecord(0, m_Buffers.size() - 1, 0);
    if (SizeInBytes != 0)
        Reserve(r, SizeInBytes);
    return r;
}

void FStorageDeviceMemoryBuf::Write(FRecord const &r, void const *pData, FStreamOffset nLength, FStreamOffset Offset) const
{
    FBuffer &Buf = GetBuf(r);
    if ( Buf.p == 0 )
        Reserve(r, Offset + nLength);

    Offset += r.BaseOffset;
    assert(Offset + nLength <= Buf.Length );
    memcpy(Buf.p + Offset, pData, nLength);
}

void FStorageDeviceMemoryBuf::Read(FRecord const &r, void *pData, FStreamOffset nLength, FStreamOffset Offset) const
{
    FBuffer &Buf = GetBuf(r);
    Offset += r.BaseOffset;
    assert(Offset + nLength <= Buf.Length );
    memcpy(pData, Buf.p + Offset, nLength);
}

void FStorageDeviceMemoryBuf::Reserve(FRecord const &r, FStreamOffset nLength) const
{
    FBuffer &Buf = GetBuf(r);
    delete Buf.p;
    Buf.p = static_cast<char*>(malloc(nLength));
    if (Buf.p == 0)
       throw std::runtime_error(fmt::format("FStorageDeviceMemoryBuf::Reserve: failed to allocate {:.2f} MB of memory", (double(nLength)/double(size_t(1) << size_t(20)))));
    Buf.Length = nLength;
}

void FStorageDeviceMemoryBuf::Delete(FRecord const &r)
{
    FBuffer &Buf = GetBuf(r);
    delete Buf.p;
    Buf.p = 0;
    Buf.Length = 0;
}


FStorageDeviceMemoryBuf::~FStorageDeviceMemoryBuf()
{
    // semi-safe...
    for ( uint i = m_Buffers.size(); i != 0; -- i )
        delete m_Buffers[i - 1].p;
    m_Buffers.clear();
}








template<class FScalar>
void FileAndMemOp(FFileAndMemOp Op, FScalar &F, FScalar *pMemBuf,
        std::size_t nMemLength, FStorageBlock const &Rec,
        FMemoryStack &Temp)
{
    // perform simple memory block/file block operations with
    // constant memory buffer size.
    std::size_t
        nSize,
        nOff = 0,
        nFileStep = std::min(size_t(Temp.MemoryLeft()), (size_t)1000000); // ! ~ 1 mb
    FScalar
        fVal = 0.0,
        *pFileBuf;
    const uint
        Byte = sizeof(FScalar);
    nFileStep = nFileStep/Byte - 1;
    assert_rt(nFileStep != 0);
    Temp.Alloc(pFileBuf, nFileStep);

    while (nOff < nMemLength) {
        nSize = std::min(nFileStep, nMemLength - nOff);
        Rec.pDevice->Read(Rec.Record, pFileBuf, nSize*Byte, nOff*Byte);

        switch( Op ){
            case OP_AddToMem:
                // add f*file_content to BufM.
                Add(pMemBuf + nOff, pFileBuf, F, nSize);
                break;
            case OP_AddToFile:
                // add f*mem_content to file
                Add(pFileBuf, pMemBuf + nOff, F, nSize);
                Rec.pDevice->Write(Rec.Record, pFileBuf, nSize*Byte, nOff*Byte);
                break;
            case OP_Dot:
                // form dot(mem,file), store result in F.
                fVal += Dot(pFileBuf, pMemBuf + nOff, nSize);
                F = fVal;
                break;
            default:
                assert(0);
        }
        nOff += nSize;
    }

    Temp.Free(pFileBuf);
}


template void FileAndMemOp(FFileAndMemOp, double&, double*, std::size_t, FStorageBlock const &, FMemoryStack &);
// template void FileAndMemOp( FFileAndMemOp, float&, float*, std::size_t, FStorageBlock const &, FMemoryStack & );

} // namespace ct
