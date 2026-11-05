// Accessibility check of the report pages (make a11y, docs/benchmarks.md).
// Usage: node benches/a11y.mjs NODE_MODULES_DIR PAGE.html...
// Each page runs axe-core in Chromium in light, dark and print (dark system
// theme, print media) at 320, 360 and 1280 px. At 320 and 360 px the
// document must not scroll sideways, with every details element closed and
// then open. Tab must reach every summary, link and scroll box outside a
// closed details element. Exit 1 on any miss.
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { join, resolve } from 'node:path';

const [modules, ...pages] = process.argv.slice(2);
const require = createRequire(join(resolve(modules), 'node_modules') + '/');
const { chromium } = require('playwright');
const axeSource = readFileSync(require.resolve('axe-core/axe.min.js'), 'utf8');
const name = (file) => file.split('/').slice(-2).join('/');

const browser = await chromium.launch();
let failures = 0;
for (const scheme of ['light', 'dark', 'print']) {
  for (const width of [320, 360, 1280]) {
    const context = await browser.newContext({
      viewport: { width, height: 900 },
      colorScheme: scheme === 'light' ? 'light' : 'dark',
      bypassCSP: true,
    });
    const page = await context.newPage();
    if (scheme === 'print') await page.emulateMedia({ media: 'print' });
    for (const file of pages) {
      await page.goto('file://' + resolve(file));
      const scrollWidth = () => page.evaluate(() => document.documentElement.scrollWidth);
      const closed = await scrollWidth();
      await page.evaluate(() => document.querySelectorAll('details').forEach((d) => (d.open = true)));
      const open = await scrollWidth();
      await page.addScriptTag({ content: axeSource });
      const axe = await page.evaluate(async () => {
        const result = await window.axe.run(document, { resultTypes: ['violations'] });
        return {
          version: window.axe.version,
          violations: result.violations.map((v) => `${v.id} (${v.nodes.length})`),
        };
      });
      const reflows = width > 360 || (closed <= width && open <= width);
      if (!reflows || axe.violations.length) failures++;
      console.log(
        `${scheme} ${width}px ${name(file)}: width closed ${closed}, open ${open};` +
          ` axe ${axe.version} violations ${axe.violations.length} ${axe.violations.join(' ')}`,
      );
    }
    await context.close();
  }
}

const context = await browser.newContext({ viewport: { width: 1280, height: 900 } });
const page = await context.newPage();
for (const file of pages) {
  await page.goto('file://' + resolve(file));
  const wanted = await page.evaluate(
    () =>
      Array.from(document.querySelectorAll('summary, .tbl[tabindex], a[href]')).filter((e) => {
        const details = e.parentElement.closest('details');
        return !details || details.open || e.tagName === 'SUMMARY';
      }).length,
  );
  const reached = new Set();
  for (let i = 0; i < 400; i++) {
    await page.keyboard.press('Tab');
    const id = await page.evaluate(() => {
      const e = document.activeElement;
      if (!e || e === document.body) return null;
      return `${e.tagName}:${Array.from(e.parentNode.children).indexOf(e)}:${e.textContent.slice(0, 40)}`;
    });
    if (id === null && reached.size) break;
    if (id) reached.add(id);
  }
  if (reached.size < wanted) failures++;
  console.log(`keyboard ${name(file)}: ${reached.size} of ${wanted} reached by Tab`);
}
await browser.close();
console.log(failures ? `a11y: ${failures} failures` : 'a11y: ok');
process.exit(failures ? 1 : 0);
