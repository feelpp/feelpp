Development FAQ
===============

### How to print data.

```cpp
auto val1 = integrate(_range=elements(mesh),_expr=print(idv(myScalar),"myScalar") );
auto val2 = integrate(_range=markedfaces(mesh,"markerName"),_expr=print(idv(myScalar),"myScalar") );
auto val2 = integrate(_range=markedfaces(mesh,"markerName"),_expr=print(trans(idv(myVector)),"myVectorTrans")*print(idv(myVector),"myVector") );
```
