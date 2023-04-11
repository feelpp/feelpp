# Offline step

```bash
./feelpp_mor_eye2brainapp --config-file /home/saigre/Documents/code/feelpp-dev/mor/examples/eye2brain/eye2brain/eye2brain.cfg
```


# Online step

- Run online
```bash
feelpp_mor_onlinerun --crbmodel.name eye2brain_p1g1
```

- Comparer RB and PFEM with Python :
```bash
/bin/python3 /home/saigre/Documents/code/feelpp-dev/mor/examples/eye2brain/compFE-RB.py -N 3
```