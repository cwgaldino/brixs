git add .
git commit -m "Documentation"

cd docs
cd build/html
tar czf /tmp/html.tgz .
cd ../../../
git checkout gh-pages
git rm -rf .
rm -r brixs
rm -r brixs.egg-info
rm -r build
rm -r dist
rm -r docs
tar xzf /tmp/html.tgz
git add .
git commit -a -m "publish the docs"
git push origin gh-pages
git checkout main
