const { chromium } = require('playwright');
(async () => {
  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage();
  
  page.on('console', msg => console.log(`BROWSER CONSOLE [${msg.type()}]: ${msg.text()}`));
  page.on('pageerror', error => console.log(`PAGE ERROR: ${error.message}`));
  
  await page.goto('http://127.0.0.1:5000');
  
  const fileInput = await page.$('#file-input');
  await fileInput.setInputFiles('/home/jpmslima/coding/Sauron/1afw.cif');
  
  await page.click('button[type="submit"]');
  console.log("Clicked submit, waiting...");
  
  await page.waitForTimeout(10000); // wait for processing
  console.log("Done waiting.");
  
  await browser.close();
})();
