kubectl exec deployment/toolbox -it -- touch /opt/chipster/toolbox/.reload/touch-me-to-reload-tools ; timeout 2 kubectl logs deployment/toolbox -f
