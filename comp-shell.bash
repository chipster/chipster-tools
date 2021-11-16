#!/bin/bash

if ! kubectl get pod | grep install-r | grep Running > /dev/null 2>&1; then
  cat <<EOF | kubectl apply -f -
apiVersion: apps/v1
kind: Deployment
metadata:
  name: install-r
spec:
  replicas: 1
  selector:
    matchLabels:
      app: install-r
  template:
    metadata:
      labels:
        app: install-r
    spec:
      containers:
      - name: install-r
        image: docker-registry.rahti.csc.fi/chipster-images-beta/comp-20.04-r-deps
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
    if kubectl get pod 2> /dev/null | grep install-r | grep Running > /dev/null; then     
      echo ""
      break
    fi
    printf "."
    sleep 1
  done
fi

echo "comp shell, press Ctrl+D to exit"
kubectl exec -it $(kubectl get pod | grep install-r | grep Running | cut -d " " -f 1) -- bash
