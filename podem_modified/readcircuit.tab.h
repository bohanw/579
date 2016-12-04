/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     GINPUT = 258,
     GOUTPUT = 259,
     GDFF = 260,
     GNOT = 261,
     GBUF = 262,
     GAND = 263,
     GOR = 264,
     GNAND = 265,
     GNOR = 266,
     EOLINE = 267,
     LPAR = 268,
     RPAR = 269,
     COMMA = 270,
     EQUAL = 271,
     NAMESTRING = 272
   };
#endif
/* Tokens.  */
#define GINPUT 258
#define GOUTPUT 259
#define GDFF 260
#define GNOT 261
#define GBUF 262
#define GAND 263
#define GOR 264
#define GNAND 265
#define GNOR 266
#define EOLINE 267
#define LPAR 268
#define RPAR 269
#define COMMA 270
#define EQUAL 271
#define NAMESTRING 272




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 35 "readcircuit.y"
{
long num;
char str[256];
GATEFUNC func;
}
/* Line 1529 of yacc.c.  */
#line 89 "readcircuit.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

