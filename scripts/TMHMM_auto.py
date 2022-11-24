from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium import webdriver
import time
import datetime

start = datetime.datetime.now()

options = webdriver.ChromeOptions()
options.binary_location = "C:\Program Files\Google\Chrome\Application\chrome.exe"
chrome_driver_binary = "C:\Program Files\Google\Chrome\Application\chromedriver.exe"
driver = webdriver.Chrome(chrome_driver_binary)
driver.get("https://services.healthtech.dtu.dk/service.php?TMHMM-2.0")
time.sleep(5)

# upload files
upload = driver.find_element_by_css_selector("input[type='file']")
upload.send_keys("D:\data\symbio\DOMAIN\TMD\actinidia_chinensis.tsv")
time.sleep(1)

# select output format
output = driver.find_element_by_css_selector("input[value='-short']").click()
time.sleep(1)

# submit
submit = driver.find_element_by_css_selector("input[type='submit']").click()
time.sleep(1)

driver.current_url
    #while True:
     #   try:
      #      driver.find_element_by_css_selector("span[title='Download the following table in tab-delimited format']").click()
       #     driver.current_url
        #    print(name+"_21.fa 分析完成")
         #   break
     #   except:
      #      now=datetime.datetime.now()
       #     cost=(now-start).seconds
        #    print("分析中，请稍等.....已花费{}秒".format(cost))
         #   time.sleep(20)

# driver.close()