#ifndef NAILTEST_H
#define NAILTEST_H

#include <stdio.h>
#include "cmplx.h"

void NailTestPrintComplexArray(const cmplx* array, int length)
{
	for (int c = 0; c < length; c++)
	{
		unsigned int lowordR = *(((int*)&array[c].r) + 0);
		unsigned int hiwordR = *(((int*)&array[c].r) + 1);
		unsigned int lowordI = *(((int*)&array[c].i) + 0);
		unsigned int hiwordI = *(((int*)&array[c].i) + 1);
		printf("0x%08x%08xr|0x%08x%08xi\n", hiwordR, lowordR, hiwordI, lowordI);
	}
}

void NailTestPrintDoubleArray(const double* array, int length)
{
	int col = 0;
	for (int c = 0; c < length; c++)
	{
		unsigned int loword = *(((int*)&array[c]) + 0);
		unsigned int hiword = *(((int*)&array[c]) + 1);
		printf("0x%08x%08x", hiword, loword);
		if (c != (length - 1))
		{
			printf(", ");
		}
		if (++col > 4)
		{
			printf("\n");
			col = 0;
		}
	}
	printf("\n");
}

void NailTestPrintFloatArray(const float* array, int length)
{
	int col = 0;
	for (int c = 0; c < length; c++)
	{
		unsigned int testhex = *((unsigned int*)&array[c]);
		printf("0x%xU", testhex);
		if (c != (length - 1))
		{
			printf(", ");
		}
		if (++col > 12)
		{
			printf("\n");
			col = 0;
		}
	}
}

void NailTestPrintIntArray(const int* array, int length)
{
	printf("new int[] { ");
	int col = 0;
	for (int c = 0; c < length; c++)
	{
		printf("%d", array[c]);
		if (c != (length - 1))
		{
			printf(",");
		}
		if (++col > 12)
		{
			printf("\n");
			col = 0;
		}
	}

	printf("}");
}

void NailTestPrintUintArray(const unsigned int* array, int length)
{
	printf("new uint[] { ");
	int col = 0;
	for (int c = 0; c < length; c++)
	{
		printf("%dU", array[c]);
		if (c != (length - 1))
		{
			printf(",");
		}
		if (++col > 12)
		{
			printf("\n");
			col = 0;
		}
	}

	printf("}");
}

void NailTestPrintShortArray(const short* array, int length)
{
	printf("new short[] { ");
	int col = 0;
	for (int c = 0; c < length; c++)
	{
		printf("%d", array[c]);
		if (c != (length - 1))
		{
			printf(",");
		}
		if (++col > 16)
		{
			printf("\n");
			col = 0;
		}
	}

	printf("}");
}

void NailTestPrintShortArrayAsInt(const short* array, int length)
{
	printf("new int[] { ");
	int col = 0;
	for (int c = 0; c < length; c++)
	{
		printf("%d", array[c]);
		if (c != (length - 1))
		{
			printf(",");
		}
		if (++col > 16)
		{
			printf("\n");
			col = 0;
		}
	}

	printf("}");
}

void NailTestPrintByteArray(const unsigned char* array, int length)
{
	printf("new byte[] { ");
	int col = 0;
	for (int c = 0; c < length; c++)
	{
		printf("%d", array[c]);
		if (c != (length - 1))
		{
			printf(",");
		}
		if (++col > 32)
		{
			printf("\n");
			col = 0;
		}
	}

	printf("}");
}

void NailTestPrintSbyteArray(const signed char* array, int length)
{
	printf("new sbyte[] { ");
	int col = 0;
	for (int c = 0; c < length; c++)
	{
		printf("%d", array[c]);
		if (c != (length - 1))
		{
			printf(",");
		}
		if (++col > 32)
		{
			printf("\n");
			col = 0;
		}
	}

	printf("}");
}

void NailTestPrintInt(char* varName, const int var)
{
	printf("%s = %d", varName, var);
}

void NailTestPrintUint(char* varName, const unsigned int var)
{
	printf("%s = 0x%xU", varName, var);
}

void NailTestPrintShort(char* varName, const short var)
{
	printf("%s = %d", varName, var);
}

void NailTestPrintSbyte(char* varName, const int var)
{
	printf("%s = %d", varName, var);
}

void NailTestPrintFloat(char* varName, const float var)
{
	printf("%s = BitConverter.ToSingle(BitConverter.GetBytes((uint)0x%xU), 0)", varName, *((unsigned int*)&var));
}


static int TestNumCounter = 0;

void NailTestPrintTestHeader(char* methodName)
{
	printf("[TestMethod]\npublic void Test_%s_%d()\n{\n", methodName, TestNumCounter++);
}

void NailTestPrintTestFooter()
{
	printf("}\n\n");
}

void NailTestPrintInputFloatArrayDeclaration(char* varName, const float* array, const int length)
{
	printf("Pointer<float> in_%s = Helpers.WrapWithArrayPointer<float>(\n", varName);
	NailTestPrintFloatArray(array, length);
	printf(");\n");
}

void NailTestPrintInputIntArrayDeclaration(char* varName, const int* array, const int length)
{
	printf("Pointer<int> in_%s = Helpers.WrapWithArrayPointer<int>(\n", varName);
	NailTestPrintIntArray(array, length);
	printf(");\n");
}

void NailTestPrintInputShortArrayDeclaration(char* varName, const short* array, const int length)
{
	printf("Pointer<short> in_%s = Helpers.WrapWithArrayPointer<short>(\n", varName);
	NailTestPrintShortArray(array, length);
	printf(");\n");
}

void NailTestPrintInputSbyteArrayDeclaration(char* varName, const signed char* array, const int length)
{
	printf("Pointer<sbyte> in_%s = Helpers.WrapWithArrayPointer<sbyte>(\n", varName);
	NailTestPrintSbyteArray(array, length);
	printf(");\n");
}

void NailTestPrintInputByteArrayDeclaration(char* varName, const unsigned char* array, const int length)
{
	printf("Pointer<byte> in_%s = Helpers.WrapWithArrayPointer<byte>(\n", varName);
	NailTestPrintByteArray(array, length);
	printf(");\n");
}

void NailTestPrintInputIntDeclaration(char* varName, const int var)
{
	printf("int in_");
	NailTestPrintInt(varName, var);
	printf(";\n");
}

void NailTestPrintInputUintDeclaration(char* varName, const unsigned int var)
{
	printf("uint in_");
	NailTestPrintUint(varName, var);
	printf(";\n");
}

void NailTestPrintInputSbyteDeclaration(char* varName, const signed char var)
{
	printf("sbyte in_");
	NailTestPrintSbyte(varName, var);
	printf(";\n");
}

void NailTestPrintInputFloatDeclaration(char* varName, const float var)
{
	printf("float in_");
	NailTestPrintFloat(varName, var);
	printf(";\n");
}

void NailTestPrintOutputFloatArrayDeclaration(char* varName, const float* array, const int length)
{
	printf("float[] expected_%s = \n", varName);
	NailTestPrintFloatArray(array, length);
	printf(";\n");
}

void NailTestPrintOutputIntArrayDeclaration(char* varName, const int* array, const int length)
{
	printf("int[] expected_%s = \n", varName);
	NailTestPrintIntArray(array, length);
	printf(";\n");
}

void NailTestPrintOutputSbyteArrayDeclaration(char* varName, const signed char* array, const int length)
{
	printf("sbyte[] expected_%s = \n", varName);
	NailTestPrintSbyteArray(array, length);
	printf(";\n");
}

void NailTestPrintOutputByteArrayDeclaration(char* varName, const unsigned char* array, const int length)
{
	printf("byte[] expected_%s = \n", varName);
	if (!array)
	{
		printf("new byte[0];\n");
		return;
	}
	NailTestPrintByteArray(array, length);
	printf(";\n");
}

void NailTestPrintOutputIntDeclaration(char* varName, int var)
{
	printf("int expected_");
	NailTestPrintInt(varName, var);
	printf(";\n");
}

void NailTestPrintOutputUintDeclaration(char* varName, unsigned int var)
{
	printf("uint expected_");
	NailTestPrintUint(varName, var);
	printf(";\n");
}

void NailTestPrintOutputShortArrayDeclaration(char* varName, const short* array, const int length)
{
	printf("short[] expected_%s = \n", varName);
	NailTestPrintShortArray(array, length);
	printf(";\n");
}

void NailTestPrintOutputShortDeclaration(char* varName, const short var)
{
	printf("short expected_");
	NailTestPrintShort(varName, var);
	printf(";\n");
}

void NailTestPrintOutputByteDeclaration(char* varName, const unsigned char var)
{
	printf("byte expected_");
	NailTestPrintInt(varName, var);
	printf(";\n");
}

void NailTestPrintOutputSbyteDeclaration(char* varName, const signed char var)
{
	printf("sbyte expected_");
	NailTestPrintInt(varName, var);
	printf(";\n");
}

void NailTestPrintOutputFloatDeclaration(char* varName, const float var)
{
	printf("float expected_");
	NailTestPrintFloat(varName, var);
	printf(";\n");
}

void NailTestPrintMemberVarInt(char* structName, char* varName, const int value)
{
	printf("%s.", structName);
	NailTestPrintInt(varName, value);
	printf(";\n");
}

void NailTestPrintMemberVarUint(char* structName, char* varName, const unsigned int value)
{
	printf("%s.", structName);
	NailTestPrintUint(varName, value);
	printf(";\n");
}

void NailTestPrintMemberVarFloat(char* structName, char* varName, const float value)
{
	printf("%s.", structName);
	NailTestPrintFloat(varName, value);
	printf(";\n");
}

void NailTestPrintMemberVarShort(char* structName, char* varName, const short value)
{
	printf("%s.", structName);
	NailTestPrintShort(varName, value);
	printf(";\n");
}

void NailTestPrintMemberVarSbyte(char* structName, char* varName, const signed char value)
{
	printf("%s.", structName);
	NailTestPrintSbyte(varName, value);
	printf(";\n");
}

void NailTestPrintMemberVarIntArray(char* structName, char* varName, const int* value, int length)
{
	if (!value)
	{
		printf("%s.%s = null;\n", structName, varName);
		return;
	}
	printf("%s.%s = new Pointer<int>(", structName, varName);
	NailTestPrintIntArray(value, length);
	printf(");\n");
}

void NailTestPrintMemberVarShortArray(char* structName, char* varName, const short* value, int length)
{
	if (!value)
	{
		printf("%s.%s = null;\n", structName, varName);
		return;
	}
	printf("%s.%s = new Pointer<short>(", structName, varName);
	NailTestPrintShortArray(value, length);
	printf(");\n");
}

void NailTestPrintMemberVarShortArrayAsInt(char* structName, char* varName, const short* value, int length)
{
	if (!value)
	{
		printf("%s.%s = null;\n", structName, varName);
		return;
	}
	printf("%s.%s = new Pointer<int>(", structName, varName);
	NailTestPrintShortArrayAsInt(value, length);
	printf(");\n");
}

void NailTestPrintMemberVarByteArray(char* structName, char* varName, const unsigned char* value, int length)
{
	if (!value)
	{
		printf("%s.%s = null;\n", structName, varName);
		return;
	}
	printf("%s.%s = new Pointer<byte>(", structName, varName);
	NailTestPrintByteArray(value, length);
	printf(");\n");
}

void NailTestPrintMemberVarSbyteArray(char* structName, char* varName, const signed char* value, int length)
{
	if (!value)
	{
		printf("%s.%s = null;\n", structName, varName);
		return;
	}
	printf("%s.%s = new Pointer<sbyte>(", structName, varName);
	NailTestPrintSbyteArray(value, length);
	printf(");\n");
}

void NailTestPrintMemberVarFloatArray(char* structName, char* varName, const float* value, int length)
{
	if (!value)
	{
		printf("%s.%s = null;\n", structName, varName);
		return;
	}
	printf("%s.%s = new Pointer<float>(", structName, varName);
	NailTestPrintFloatArray(value, length);
	printf(");\n");
}

#endif