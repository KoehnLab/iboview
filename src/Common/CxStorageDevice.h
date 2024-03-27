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

#ifndef CX_STORAGE_DEVICE_H
#define CX_STORAGE_DEVICE_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "CxTypes.h"

namespace ct {

typedef uint64_t
    FStreamSize;
typedef uint64_t
    FStreamOffset;
uint64_t const
    UnknownSize = static_cast<FStreamOffset>(-1);

struct FStorageBlock;
struct FMemoryStack;

/// An object managing persistent storage. Abstracts real or virtual file systems.
/// Intended for interface compatibility with Molpro/ITF.
struct FStorageDevice
{
    typedef FStreamOffset
        FOffset;

    // identifies storage records within the device.
    // Resource is owned by the device, not by the FRecord structure.
    struct FRecord
    {
        uint
            iFile,
            iRecord;
        FStreamOffset
            BaseOffset,
            EndOffset;
        uint
            Id; // for internal use of storage device
        FRecord() : BaseOffset(0), EndOffset(UnknownSize), Id(0) {}
        FRecord(uint iFile_, uint iRecord_ = 0, FStreamOffset nBaseOffset_ = 0, FStreamOffset nLength_ = UnknownSize )
            : iFile(iFile_), iRecord(iRecord_), BaseOffset(nBaseOffset_), EndOffset(UnknownSize), Id(0)
        {
            if ( nLength_ != UnknownSize)
                SetLength(nLength_);
        }

        void SetLength(FStreamOffset NewLength) const;
    };


    // allocate a new record for temporary data.
    virtual FRecord AllocNewRecord(FStreamOffset SizeInBytes) = 0;
    FStorageBlock AllocNewBlock(FStreamOffset SizeInBytes);

    virtual void Write( FRecord const &r, void const *pData, FStreamOffset nLength, FStreamOffset Offset = 0 ) const = 0;
    virtual void Read( FRecord const &r, void *pData, FStreamOffset nLength, FStreamOffset Offset = 0 ) const = 0;
    virtual void Reserve( FRecord const &r, FStreamOffset nLength ) const = 0;
    virtual void Delete( FRecord const &r ) = 0;

    virtual ~FStorageDevice(); // = 0;

    FStorageDevice() {};
private:
    FStorageDevice(FStorageDevice const &other); // not implemented, don't copy!
    void operator = (FStorageDevice const &other); // not implemented, don't copy!
};

typedef FStorageDevice::FRecord FRecord;

/// Convenience-combination of storage device and block. Note that these
/// objects DO NOT OWN the data they refer to. The data belongs to the storage
/// device. As such storage blocks are copy-able and don't contain state.
struct FStorageBlock{
    FStorageDevice
        *pDevice;
    FStorageDevice::FRecord
        Record;

    FStorageBlock() : pDevice(0) {};
    FStorageBlock( FStorageDevice &Device_, FStorageDevice::FRecord const &Record_ )
        : pDevice(&Device_), Record( Record_ )
    {}

    FStorageBlock( FStorageDevice const &Device_, FStorageDevice::FRecord const &Record_ )
        : pDevice(const_cast<FStorageDevice*>(&Device_)), Record( Record_ )
    {}             // ^- yeah, I know.

    void RawWrite( void const *pData, FStreamOffset nLength, FStreamOffset Offset = 0 ) {
        pDevice->Write(Record, pData, nLength, Offset);
    }
    void RawRead( void *pData, FStreamOffset nLength, FStreamOffset Offset = 0 ) const {
        pDevice->Read(Record, pData, nLength, Offset);
    }

    template<class FScalar>
    void Write( FScalar const *pData, FStreamOffset nLength, FStreamOffset Offset = 0 ) {
        RawWrite(pData, sizeof(FScalar) * nLength, sizeof(FScalar) * Offset);
    }

    template<class FScalar>
    void Read( FScalar *pData, FStreamOffset nLength, FStreamOffset Offset = 0 ) const {
        RawRead(pData, sizeof(FScalar) * nLength, sizeof(FScalar) * Offset);
    }

    void Delete() {
        if ( pDevice != 0 )
            pDevice->Delete(Record);
    }
};

std::ostream &operator << ( std::ostream &out, FStorageDevice::FRecord const &r );
std::ostream &operator << ( std::ostream &out, FStorageBlock const &r );


/// Implements FStorageDevice interface based on basic posix file system interfaces.
/// Each record is represented by a file; files are put into locations returned by ::tmpfile().
struct FStorageDevicePosixFs : public FStorageDevice
{
    // create a record for an existing file handle. Takes ownership---fclose will be called on this if the storage device goes out of scope.
    FRecord AllocNewRecord(FStreamOffset SizeInBytes, FILE* FileHandle);
    // creates a record connected to a file from tmpfile(). This may end up in a directory with little space, however.
    FRecord AllocNewRecord(FStreamOffset SizeInBytes); // override.
    FRecord AllocNewRecord_TempFileInDir(FStreamOffset SizeInBytes, std::string const &Directory, std::string const &FileNameTemplate = "storage-{}.tmp");
    void Write( FRecord const &r, void const *pData, FStreamOffset nLength, FStreamOffset Offset = 0 ) const; // override
    void Read( FRecord const &r, void *pData, FStreamOffset nLength, FStreamOffset Offset = 0 ) const; // override
    void Reserve( FRecord const &r, FStreamOffset nLength ) const; // override
    void Delete( FRecord const &r ); // override

    FStorageDevicePosixFs();
    ~FStorageDevicePosixFs();
protected:
    typedef std::map<uint, FILE*>
        FFileIdMap;
    FFileIdMap
        m_FileIds;
    typedef std::set<FILE*>
        FFileHandleSet;
    FFileHandleSet
        m_ControlledFiles; // to make sure that we don't register one multiple times...
    FILE *GetHandle(FRecord const &r) const;
};


/// Implements FStorageDevice with data kept in memory.
/// Each record is represented by a in-memory storage block.
struct FStorageDeviceMemoryBuf : public FStorageDevice
{
    FRecord AllocNewRecord(FStreamOffset SizeInBytes = 0);
    void Write( FRecord const &r, void const *pData, FStreamOffset nLength, FStreamOffset Offset = 0 ) const;
    void Read( FRecord const &r, void *pData, FStreamOffset nLength, FStreamOffset Offset = 0 ) const;
    void Reserve( FRecord const &r, FStreamOffset nLength ) const; // override
    void Delete( FRecord const &r );

    virtual ~FStorageDeviceMemoryBuf();
private:
    struct FBuffer {
        char *p;
        std::size_t Length;
    };
    typedef std::vector<FBuffer>
        FBufferSet;
    FBufferSet
        m_Buffers;
    FBuffer &GetBuf(FRecord const &r) const;
};

enum FFileAndMemOp {
    OP_AddToMem,
    OP_Dot,
    OP_AddToFile,
    OP_WriteFile
};

// provides some simple constant buffer size file+memory operations (see .cpp).
// This template is bound for 64bit and 32bit floats.
template<class FScalar>
void FileAndMemOp( FFileAndMemOp Op, FScalar &f, FScalar *pMemBuf,
    std::size_t nMemLength, FStorageBlock const &Rec,
    FMemoryStack &Temp );

} // namespace ct

#endif // CX_STORAGE_DEVICE_H
