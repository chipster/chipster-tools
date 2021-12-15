#!/bin/bash

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -gt 1 ]; then
  echo "Open shell in container image"
  echo "Usage: $0 [IMAGE]"
  exit 0
fi

image=${1:-comp-20.04-r-deps}
name="comp-shell"

if ! kubectl get pod | grep $name | grep Running > /dev/null 2>&1; then
  cat <<EOF | kubectl apply -f -
apiVersion: apps/v1
kind: Deployment
metadata:
  name: $name
spec:
  replicas: 1
  selector:
    matchLabels:
      app: $name
  template:
    metadata:
      labels:
        app: $name
    spec:
      containers:
      - name: $name
        image: docker-registry.rahti.csc.fi/chipster-images-beta/$image
        volumeMounts:
        - name: tools-bin
          mountPath: /mnt/tools
      volumes:
      - name: tools-bin
        hostPath:
          path: /mnt/tools-bin
          type: Directory
EOF

  echo "waiting pod to start"
  for i in $(seq 30); do
    if kubectl get pod 2> /dev/null | grep $name | grep Running > /dev/null; then     
      echo ""
      break
    fi
    printf "."
    sleep 1
  done
fi

echo "comp shell, press Ctrl+D to exit"
kubectl exec -it $(kubectl get pod | grep $name | grep Running | cut -d " " -f 1) -- bash

kubectl delete deployment $name
