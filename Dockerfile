FROM ubuntu:22.04 as builder

RUN apt -y update -qq && apt -y upgrade && \
	DEBIAN_FRONTEND=noninteractive apt -y install \
	autoconf automake gawk gcc g++ clang make software-properties-common vim 

WORKDIR /usr/src 

# Put the code into a subdir, so we don't copy the Makefile and 
# Dockerfile into the container, they are not needed here.
COPY . .

RUN autoupdate && autoreconf && ./configure && make && make install

FROM ubuntu:22.04 as runtime

COPY --from=builder /usr/local/bin/rfmix /usr/local/bin/simulate \
	/usr/local/bin/

# we map the user owning the image so permissions for any mapped 
# input/output paths set by the user will work correctly
CMD rfmix
