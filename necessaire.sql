


------############################################-------
------               COMMANDE GIT                 -------
------############################################-------

-- regarder le contenu d un fichier sur remote git
git show origin/main:requirements.txt



------############################################-------
------               COMMANDE STREAMLIT           -------
------############################################-------


--list process ki tourne
pgrep -f streamlit
--kill tout les process
pkill -f streamlit
--kill un process
kill -9 3058


------############################################-------
------               COMMANDE SSH key             -------
------############################################-------



-- verifier que l agent tourne 
commad = eval "$(ssh-agent -s)"
result = Agent pid 22695


--lister les clé qui existent
ls -la ~/.ssh/

--ajouter une clé ( la clé c id_ed25519 )
ssh-add ~/.ssh/id_ed25519 

--verifier si elle bien ajoutée
ssh-add -l


--tester la connexion git
ssh -T git@github.com

--puis on peut poushé