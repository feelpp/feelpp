ARG IMAGE=ubuntu:jammy
FROM ${IMAGE}

ARG FLAVOR=ubuntu
ARG DIST=jammy
ARG VERSION=0.110.0
ARG CHANNEL
ENV CHANNEL=${CHANNEL:-latest}


RUN if test "$FLAVOR" = "ubuntu"; then \
        sed -i 's/http:\/\/archive.ubuntu.com\/ubuntu\//http:\/\/miroir.univ-lorraine.fr\/ubuntu\//' /etc/apt/sources.list; \
    fi
RUN if test "$FLAVOR" = "debian"; then \
        NEW_MIRROR="http://miroir.univ-lorraine.fr/debian/" ;\
        sed -i "s|URIs: .*|URIs: $NEW_MIRROR|" /etc/apt/sources.list.d/debian.sources;\
        fi
RUN apt-get update \
        && apt-get install -y \
        vim \
        ssh \
        sudo \
        wget gpg git \
        software-properties-common ;\
        rm -rf /var/lib/apt/lists/*

RUN useradd --user-group --create-home --shell /bin/bash feelpp ;\
        echo "feelpp ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

ENV DEBIAN_FRONTEND=noninteractive        
RUN if test "$FLAVOR" = "ubuntu"; then \
        echo "deb http://miroir.univ-lorraine.fr/ubuntu $DIST-proposed main restricted" |  tee -a /etc/apt/sources.list; \ 
        echo "deb http://miroir.univ-lorraine.fr/ubuntu $DIST-proposed universe" |  tee -a /etc/apt/sources.list; \
        echo "deb http://miroir.univ-lorraine.fr/ubuntu $DIST-proposed multiverse" |  tee -a /etc/apt/sources.list; \
    fi; \
    apt-get -qq update; \
    apt-get install -y tzdata dos2unix;\
    rm -rf /var/lib/apt/lists/*

RUN wget -qO - http://apt.feelpp.org/apt.gpg | gpg --dearmor > feelpp.gpg; \
        install -o root -g root -m 644 feelpp.gpg /etc/apt/trusted.gpg.d; \
        echo "deb [signed-by=/etc/apt/trusted.gpg.d/feelpp.gpg] http://apt.feelpp.org/$FLAVOR/$DIST $DIST ${CHANNEL}" | tee -a /etc/apt/sources.list.d/feelpp.list; \
        rm -f feelpp.gpg; \
        apt-get update ;\
        packages=$(apt-cache search '^feelpp-*' | grep -v 'dbgsym' | awk '{print $1}');\
        apt-get install -y --no-install-recommends \
                         g++ clang man-db make cmake rsync ipython3 \
                         python3-pip \
                         $packages 'libfeelpp1' 'python3-feelpp*'; \
        rm -rf /var/lib/apt/lists/* ;

# install venv
RUN apt-get update; \
    apt-get install -y python3-venv; \
    rm -rf /var/lib/apt/lists/* ;

RUN     echo "export FEELPP_REPOSITORY=/feel" >> ~feelpp/.bashrc;\
        echo "export OMPI_MCA_btl_vader_single_copy_mechanism=none" >> ~feelpp/.bashrc

USER feelpp


RUN    python3 -m venv /home/feelpp/.venv \
   &&  . /home/feelpp/.venv/bin/activate \
   &&  pip3 install --no-cache-dir --upgrade pip setuptools wheel pyvista xvfbwrapper ipykernel plotly pandas trame trame-code tabulate

CMD ["/bin/bash", "-c", "source /home/feelpp/.venv/bin/activate && exec /bin/bash"]
