CC=gcc
# to build on sundance: CC=gcc -mcpu=v9 -m64
ifeq (${COPT},)
    COPT=-O -g
endif
CFLAGS=
HG_DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE}

#global external libraries 
L=

#default to not using ssl
#ifeq (${USE_SSL},)
    USE_SSL=0
#endif

#ifeq (${USE_SSL},1)
#    L+=-lssl
#    HG_DEFS+=-DUSE_SSL
#endif

#HG_INC=-I../inc -I../../inc -I../../../inc -I../../../../inc -I../../../../../inc
HG_INC=-Iinc -I../inc

#default to not using bam
ifeq (${USE_BAM},)
    USE_BAM=0
endif

ifeq (${USE_BAM},1)
    ifeq (${SAMPATH},)
      ifeq (${MACHTYPE},x86_64)
        SAMPATH = /hive/data/outside/samtools/samtools
      else
        SAMPATH = /hive/data/outside/samtools/samtools/${MACHTYPE}
      endif
    endif
    HG_INC += -I${SAMPATH}
    HG_DEFS+=-DUSE_BAM
endif

ifeq (${HG_WARN},)
  ifeq (darwin,$(findstring darwin,${OSTYPE}))
      HG_WARN = -Wall -Wno-unused-variable -Wno-long-double
      HG_WARN_UNINIT=
  else
    ifeq (solaris,$(findstring solaris,${OSTYPE}))
      HG_WARN = -Wall -Wformat -Wimplicit -Wreturn-type
      HG_WARN_UNINIT=-Wuninitialized
    else
      HG_WARN = -Wall -Werror -Wformat -Wimplicit -Wreturn-type
      HG_WARN_UNINIT=-Wuninitialized
    endif
  endif
  # -Wuninitialized generates a warning without optimization
  ifeq ($(findstring -O,${COPT}),-O)
     HG_WARN += ${HG_WARN_UNINIT}
  endif
endif

# this is to hack around many make files not including HG_WARN in
# the link line
CFLAGS += ${HG_WARN}

ifeq (${SCRIPTS},)
    SCRIPTS=/cluster/bin/scripts
endif
ifeq (${CGI_BIN},)
    CGI_BIN=/usr/local/apache/cgi-bin
endif
ifeq (${DOCUMENTROOT},)
    DOCUMENTROOT=/usr/local/apache/htdocs
endif
ifeq (${BINDIR},)
    BINDIR = ${HOME}/bin/${MACHTYPE}
endif
ifeq (${ENCODE_PIPELINE_BIN},)
    ENCODE_PIPELINE_BIN=/cluster/data/encode/pipeline/bin
endif

MKDIR=mkdir -p
ifeq (${STRIP},)
   STRIP=strip
endif
CVS=cvs

# portable naming of compiled executables: add ".exe" if compiled on 
# Windows (with cygwin).
ifeq (${OS}, Windows_NT)
  AOUT=a.exe
  EXE=.exe
else
  AOUT=a.out
  EXE=
endif

# location of stringify program
STRINGIFY = ${BINDIR}/stringify

%.o: %.c
	${CC} ${COPT} ${CFLAGS} ${HG_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -o $@ -c $<
