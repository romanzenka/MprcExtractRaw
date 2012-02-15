#pragma once

#include <io.h>
#include <iostream>

typedef float TJavaFloat;
typedef __int64 TJavaLong;
typedef __int32 TJavaInt;

// Simple output buffer that does the endianness conversion to produce values compatible with Java.
// All functions are in the header and inlined for performance
class OutputBuffer
{
private:
	char *buffer;
	int bufferSize;
	int bufferUsed;

public:
	OutputBuffer() { 
		buffer=NULL; 
		bufferSize=0; 
		bufferUsed=0; 
	}

	~OutputBuffer() {
		if(buffer!=NULL) {
			delete [] buffer;
		}
	}

	inline char* get() {
		return buffer;
	}

	inline int usedSize() {
		return bufferUsed;
	}

	inline void addJavaLong(__int64 value) {
		ensureBufferHasSpace(8);
		char *current = buffer+bufferUsed;
		char *valuePtr = (char*)&value;
		// Swap the endianness
		*current++ = valuePtr[7];
		*current++ = valuePtr[6];
		*current++ = valuePtr[5];
		*current++ = valuePtr[4];
		*current++ = valuePtr[3];
		*current++ = valuePtr[2];
		*current++ = valuePtr[1];
		*current++ = valuePtr[0];
		bufferUsed+=8;
	}
	inline void addJavaInt(__int32 value) {
		ensureBufferHasSpace(4);
		char *current = buffer+bufferUsed;
		char *valuePtr = (char*)&value;
		// Swap the endianness
		*current++ = valuePtr[3];
		*current++ = valuePtr[2];
		*current++ = valuePtr[1];
		*current++ = valuePtr[0];
		bufferUsed+=4;
	}
	inline void addJavaFloat(float value) {
		ensureBufferHasSpace(4);
		char *current = buffer+bufferUsed;
		char *valuePtr = (char*)&value;
		// Swap the endianness
		*current++ = valuePtr[3];
		*current++ = valuePtr[2];
		*current++ = valuePtr[1];
		*current++ = valuePtr[0];
		bufferUsed+=4;
	}

	// Makes sure that the buffer would be able to fit given amount of bytes
	inline void ensureBufferHasSpace(int freeSpace) {
		if(buffer==NULL) { reallocateBuffer(freeSpace*2); return; }
		reallocateBuffer(bufferUsed+freeSpace*2);
	}

	inline void clear() {
		bufferUsed=0;
	}

private:
	inline void reallocateBuffer(int newSpace) {
		if(newSpace<=bufferSize) return;		
		char *newBuffer = new char[newSpace];
		memcpy(newBuffer, buffer, bufferSize);		
		bufferSize = newSpace;		
		delete [] buffer;
		buffer=newBuffer;
	}
};
