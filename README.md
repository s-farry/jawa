# jawa
#procedure for tagging
git checkout -b vX.X develop
git commit -a -m "new version"
git checkout master
git merge --no-ff vX.X
git tag -a X.X
git push origin master
git push origin X.X
#optional - back to develop
git checkout develop
