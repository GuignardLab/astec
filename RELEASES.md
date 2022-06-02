[Retour à l'index](./README.md)

# Étapes à suivre pour faire une release d'astec


## release astec

- vérifier que les branches locales de *master* et *develop* soient bien à jour

```
git flow release start 2.8.0
```

  - pour une release (nouvelles features/ changement d'API), le **major** et/ou le **minor** sont changés avec le **patch level** à 0 (zéro)
  - pour une release de bugfix/hotfix le **patch level** est incrémenté seul

- mettre à jour le **CHANGES.md** en ajoutant en haut de fichier ce que contient la release (sous la ligne **## NEXT VERSION**)

```
## version 2.8.0 - 2020-04-06
- Add setBounds methods for range parameter and fix bugs
- Register std::array to QMETATYPE system
- ...
```

- mettre à jour la version dans le **setup.py** global

```
    ...
    version="2.8.0",
    ...
```

- mettre à jour la version pour la documentation dans **doc/source/conf.py** 

```
release = '2.8.0'
```

- commit local des mises à jour avec le commentaire `release x.y.z`
- une fois la release terminée (```git flow release finish x.y.z```):
- ```git push --follow-tags origin develop master``` (normalement, cela pousse les branches develop et master, et le tag)

