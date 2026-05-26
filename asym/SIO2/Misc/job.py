import requests
from bs4 import BeautifulSoup

file= "https://docs.google.com/document/d/e/2PACX-1vSvM5gDlNvt7npYHhp_XfsJvuntUhq184By5xO_pA4b_gCWeXb6dM6ZxwN8rE6S4ghUsCj2VKR21oEP/pub"


def print_message(doc_url: str) ->None:
    response = requests.get(doc_url)
    response.raise_for_status()
    
    soup= BeautifulSoup(response.text, "html.parser")
    points=[]
  
    for row in soup.find_all("tr"):
      cells = [cell.get_text(strip=True) for cell in row.find_all(["td", "th"])]
      if len(cells) < 3:
          continue
  
      try:
          x= int(cells[0])
          char = cells[1] if cells[1] else " "
          y= int(cells[2])
          points.append((x, y, char))
      except ValueError:
        continue
        
    if not points:
      raise ValueError("No coordinate data")
  
    max_x = max(x for x, _, _ in points)
    max_y = max(y for _, y, _ in points)

    grid = [[" " for _ in range(max_x + 1)] for _ in range(max_y + 1)]
    
    for x, y, char in points:
        grid[y][x]= char
    
    for y in range(max_y, -1, -1):
        print("".join(grid[y]).rstrip())

print_message(file)