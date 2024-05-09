Suite du Rendu

On a pu obtenir une meillieure performance avec 12 threads avec la taille du blocking : 2x8x100

L'obetention de ces valeurs est grace à l'etude de "stencil shape", elle se differencie de notre 1ere taille de blocking (8x8x8) qui est basée sur "caches size".

Pour avoir une idee brieve (schemas ne sont pas exactes) : 

En effet, pour 1 iteration de la boule externe (ie: x fixe) de la version ancienne (xyz->Jacobi):

<img width="508" alt="xytjacob" src="https://github.com/yutaizhao/TOP-project/assets/15853429/897d5886-5304-4fe6-bab3-36f777343949">

Or, pour 1 iteration de la boule externe (ie: jacobi fixe) de la nouvelle version (Jacobi->xyz):

<img width="427" alt="jacobxyz" src="https://github.com/yutaizhao/TOP-project/assets/15853429/39ee131e-4b05-4be3-9049-0f2203ddc17d">

Ainsi, pour la 1ere version, on essaye d'avoir tous les elements de xyz en gathering les cellules voisins ...mais 
pour la 2eme version, on essaye de beneficie au maximum le parcours sur les elements de l'axe z puis on ajuste y puis x , car le parcours de la nouvelle version forme des plans en profondeur(l'axe z) 

On a enfin :

![final latency](https://github.com/yutaizhao/TOP-project/assets/15853429/4eb3dc81-17c9-477c-bca1-f6a306f67b3c)

![final scalability](https://github.com/yutaizhao/TOP-project/assets/15853429/bc64d35e-7eb2-4e0a-a2bc-e2dc8d2d5ba5)

