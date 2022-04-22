# 1
from PyPDF2 import PdfFileReader
pdfReader = PdfFileReader('E:/HPK1/HPK1_patent/29.WO2016_090300.pdf', 'rb')
first_page = pdfReader.getPage(5)
print(first_page.extractText())

for pg in range(0, pdfReader.numPages):
    page = pdfReader.getPage(pg)
    print(page.extractText())

# 2
from PyPDF2 import PdfFileReader
pdfReader = PdfFileReader('E:/발표자료/Hgraph/data/E_Combichem/Concepts of Artificial Intelligence for Computer-Assisted Drug.pdf', 'rb')
first_page = pdfReader.getPage(5)
print(first_page.extractText())

for pg in range(0, pdfReader.numPages):
    page = pdfReader.getPage(pg)
    print(page.extractText())

# 3
from tika import parser
raw = parser.from_file('E:/HPK1/HPK1_patent/29.WO2016_090300.pdf')
print(raw['content'])

# 4
from pdfminer3.layout import LAParams, LTTextBox
from pdfminer3.pdfpage import PDFPage
from pdfminer3.pdfinterp import PDFResourceManager
from pdfminer3.pdfinterp import PDFPageInterpreter
from pdfminer3.converter import PDFPageAggregator
from pdfminer3.converter import TextConverter
import io

resource_manager = PDFResourceManager()
fake_file_handle = io.StringIO()
converter = TextConverter(resource_manager, fake_file_handle, laparams=LAParams())
page_interpreter = PDFPageInterpreter(resource_manager, converter)

with open('E:/HPK1/HPK1_patent/29.WO2016_090300.pdf', 'rb') as fh:

    for page in PDFPage.get_pages(fh,
                                  caching=True,
                                  check_extractable=True):
        page_interpreter.process_page(page)

    text = fake_file_handle.getvalue()

# close open handles
converter.close()
fake_file_handle.close()
print(text)

