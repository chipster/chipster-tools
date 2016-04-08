FROM java:8

RUN mkdir /opt/chipster-tools
RUN chown -R chipster:chipster /opt/chipster-*

USER 1001

# no need to build these projects, clone directly to /opt
RUN cd /opt; git clone --branch ${TOOLS_BRANCH:-openshift} https://github.com/chipster/chipster-tools.git --depth=1

WORKDIR /opt/chipster-web-server  

CMD java -cp lib/*: ${JAVA_CLASS}

